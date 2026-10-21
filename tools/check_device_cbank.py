#!/usr/bin/env python3
"""Check a device-linked library's compiler-constant references.

    check_device_cbank.py [--cuobjdump PATH] [--remove-on-failure] LIB OBJ [OBJ ...]

LIB is the linked shared library (libterco.so); the OBJs are the objects it
was linked from. Every kernel reads the compiler's constants -- double
literals, and the jump tables of `select case` -- from constant bank 2
(c[0x2]). Each function's constants start out in its own
`.nv.constant2.<fn>` section and nvlink lays them out, per kernel, into
that kernel's bank, patching the c[0x2] offset in every instruction that
reads one (relocations R_CUDA_CONST_FIELD19_40 for c[0x2][imm] and
R_CUDA_CONST_FIELD21_38 for the jump tables' c[0x2][Rx+imm]).

A bank holds 64 KiB. When the constants a kernel can reach add up to more
than that, nvlink (CUDA 12.9, the last toolkit for sm_70) de-duplicates the
constants of every kernel in the image into one shared bank, and in that
path it leaves some jump-table offsets unpatched: the kernel then loads its
branch targets from someone else's constants and either faults
(CUDA_ERROR_INVALID_PC / ILLEGAL_INSTRUCTION) or jumps to a wrong case and
gives silently wrong numbers. nvlink prints nothing. The meta-GGA dispatch
reached this first (trc_xc_xc_points needed 35808 bytes before it).

So this script fails the build when either holds:
  * a c[0x2] reference in the linked code does not point at the bytes its
    symbol had in the object (checked in the bank of every kernel that
    reaches the function), or
  * every kernel shares one byte-identical bank -- the over-64-KiB path,
    whether or not it happened to corrupt a reference this time.
It prints the largest per-kernel bank, the margin left under 64 KiB.

Needs cuobjdump from the CUDA toolkit nvfortran builds against. Standard
library only.
"""
import argparse
import collections
import concurrent.futures
import os
import shutil
import struct
import subprocess
import sys
import tempfile

SHT_RELA, SHT_REL, SHT_SYMTAB = 4, 9, 2
R_CUDA_ABS32_LO_32, R_CUDA_ABS32_HI_32, R_CUDA_ABS47_34 = 0x38, 0x39, 0x3a
R_CUDA_CONST_FIELD19_40, R_CUDA_CONST_FIELD21_38 = 0x40, 0x42
BANK_LIMIT = 65536


class Elf:
    """The few parts of an ELF64 (CUDA cubin) this check needs."""

    def __init__(self, path):
        d = self.d = open(path, 'rb').read()
        shoff, = struct.unpack_from('<Q', d, 0x28)
        shentsize, shnum, shstrndx = struct.unpack_from('<HHH', d, 0x3a)
        raw = [struct.unpack_from('<IIQQQQIIQQ', d, shoff + i * shentsize) for i in range(shnum)]
        names = raw[shstrndx][4]
        self.secs = [dict(name=self._cstr(names + r[0]), type=r[1], off=r[4], size=r[5],
                          link=r[6], info=r[7]) for r in raw]
        self.syms = []
        for s in self.secs:
            if s['type'] != SHT_SYMTAB:
                continue
            strtab = self.secs[s['link']]['off']
            for i in range(s['size'] // 24):
                nm, _, _, shndx, val, sz = struct.unpack_from('<IBBHQQ', d, s['off'] + 24 * i)
                self.syms.append(dict(name=self._cstr(strtab + nm), shndx=shndx, value=val, size=sz))

    def _cstr(self, o):
        return self.d[o:self.d.index(b'\0', o)].decode()

    def data(self, i):
        s = self.secs[i]
        return self.d[s['off']:s['off'] + s['size']]

    def relocs(self):
        """(target section, offset, type, symbol, addend) for every relocation."""
        out = []
        for s in self.secs:
            if s['type'] not in (SHT_REL, SHT_RELA):
                continue
            ent = 24 if s['type'] == SHT_RELA else 16
            for i in range(s['size'] // ent):
                if ent == 24:
                    off, info, add = struct.unpack_from('<QQq', self.d, s['off'] + ent * i)
                else:
                    (off, info), add = struct.unpack_from('<QQ', self.d, s['off'] + ent * i), 0
                out.append((s['info'], off, info & 0xffffffff, self.syms[info >> 32], add))
        return out


def bank_offset(ins, rtype, addend):
    """The c[0x2] byte offset a linked 128-bit instruction reads."""
    if rtype == R_CUDA_CONST_FIELD21_38:  # 5-bit bank, signed 16-bit byte offset
        imm = (ins >> 38) & 0xffff
        return imm - ((imm & 0x8000) << 1) - addend
    return ((ins >> 40) & 0x3fff) * 4 - addend  # 5-bit bank, 14-bit word offset


def extract(cuobjdump, path, dest):
    """Every device ELF in `path`, as {arch: cubin path}."""
    os.makedirs(dest, exist_ok=True)
    subprocess.run([cuobjdump, '-xelf', 'all', os.path.abspath(path)], cwd=dest,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    found = {}
    for f in os.listdir(dest):
        if f.endswith('.cubin'):
            arch = f.rsplit('.', 2)[-2]
            found[arch] = os.path.join(dest, f)
    return found


def check_arch(linked, objects):
    L = Elf(linked)
    text = collections.defaultdict(list)
    for i, s in enumerate(L.secs):
        if s['name'].startswith('.text.'):
            text[s['name'][6:]].append(i)
    banks = {s['name'][len('.nv.constant2.'):]: i for i, s in enumerate(L.secs)
             if s['name'].startswith('.nv.constant2.')}
    lrel = collections.defaultdict(list)
    calls = collections.defaultdict(set)
    for r in L.relocs():
        tgt, off, rt, sym, _ = r
        lrel[tgt].append(r)
        if rt in (R_CUDA_ABS47_34, R_CUDA_ABS32_LO_32, R_CUDA_ABS32_HI_32) \
                and sym['shndx'] < len(L.secs) and L.secs[sym['shndx']]['name'].startswith('.text.'):
            calls[tgt].add(sym['shndx'])
    reach = collections.defaultdict(set)
    for k in banks:
        if len(text.get(k, [])) != 1:
            continue
        stack, seen = [text[k][0]], set()
        while stack:
            x = stack.pop()
            if x not in seen:
                seen.add(x)
                reach[x].add(k)
                stack.extend(calls[x])

    sizes = {k: L.secs[i]['size'] for k, i in banks.items()}
    contents = {L.data(i) for i in banks.values()}
    merged = len(banks) > 1 and len(contents) == 1
    checked, wrong = 0, []
    for path in objects:
        O = Elf(path)
        orel = O.relocs()
        # Object-local duplicates (the CUDA math helpers every object carries a
        # copy of) are matched to their linked copy through the call sites.
        omap = {}
        for tgt, off, rt, sym, _ in orel:
            fn = O.secs[tgt]['name'][6:]
            if rt == R_CUDA_ABS47_34 and len(text.get(fn, [])) == 1:
                at = {o: s['shndx'] for _, o, t, s, _ in lrel[text[fn][0]] if t == R_CUDA_ABS47_34}
                if off in at:
                    omap[sym['shndx']] = at[off]
        for tgt, off, rt, sym, add in orel:
            if rt not in (R_CUDA_CONST_FIELD19_40, R_CUDA_CONST_FIELD21_38):
                continue
            cands = text.get(O.secs[tgt]['name'][6:], [])
            if len(cands) != 1:
                cands = [omap[tgt]] if tgt in omap else []
            if not cands or not reach.get(cands[0]):
                continue  # inlined everywhere or never reached from a kernel
            ins = int.from_bytes(L.data(cands[0])[off:off + 16], 'little')
            x = bank_offset(ins, rt, add)
            n = max(4, min(sym['size'] or 8, 8))
            want = O.data(sym['shndx'])[sym['value']:sym['value'] + n]
            checked += 1
            bad = [k for k in reach[cands[0]] if L.data(banks[k])[x:x + n] != want]
            if bad:
                wrong.append((O.secs[tgt]['name'][6:], off, sym['name'], len(bad), sorted(bad)[0]))
    return dict(kernels=len(banks), sizes=sizes, merged=merged, checked=checked, wrong=wrong)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--cuobjdump', default=shutil.which('cuobjdump') or 'cuobjdump')
    ap.add_argument('--remove-on-failure', action='store_true',
                    help='delete LIB when the check fails, so a broken library cannot be used')
    ap.add_argument('lib')
    ap.add_argument('objects', nargs='+')
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix='trc_cbank_')
    try:
        linked = extract(a.cuobjdump, a.lib, os.path.join(tmp, 'lib'))
        if not linked:
            print('check_device_cbank: no device code in %s, nothing to check' % a.lib)
            return 0
        objs = collections.defaultdict(list)
        with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as ex:
            futs = [ex.submit(extract, a.cuobjdump, o, os.path.join(tmp, 'o%05d' % i))
                    for i, o in enumerate(a.objects)]
            for f in futs:
                for arch, p in f.result().items():
                    objs[arch].append(p)
        failed = False
        for arch, cub in sorted(linked.items()):
            r = check_arch(cub, objs.get(arch, []))
            big = max(r['sizes'].items(), key=lambda kv: kv[1]) if r['sizes'] else ('-', 0)
            if r['merged']:
                print('check_device_cbank: %s: %d kernels, %d c[0x2] references checked'
                      % (arch, r['kernels'], r['checked']))
            else:
                print('check_device_cbank: %s: %d kernels, %d c[0x2] references checked, '
                      'largest bank %d bytes (%.0f%% of 64 KiB) in %s'
                      % (arch, r['kernels'], r['checked'], big[1], 100.0 * big[1] / BANK_LIMIT, big[0]))
            if r['merged']:
                failed = True
                print('check_device_cbank: ERROR: %s: all %d kernels share one %d-byte constant bank: '
                      'some kernel reaches more than 64 KiB of compiler constants and nvlink took '
                      'its de-duplicating path, which leaves jump-table offsets unpatched. Shrink '
                      'what that kernel reaches (split its dispatch, or compile the heavy sources '
                      'with -Wptxas,--disable-optimizer-constants).' % (arch, r['kernels'], big[1]))
            for fn, off, sym, nbad, k in r['wrong'][:20]:
                print('check_device_cbank: ERROR: %s: %s+0x%x reads the wrong c[0x2] bytes for %s '
                      '(in %d kernel(s), e.g. %s)' % (arch, fn, off, sym, nbad, k))
            if r['wrong']:
                failed = True
                print('check_device_cbank: ERROR: %s: %d constant references point at the wrong data'
                      % (arch, len(r['wrong'])))
        if failed:
            if a.remove_on_failure and os.path.exists(a.lib):
                os.remove(a.lib)
                print('check_device_cbank: removed %s' % a.lib)
            return 1
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())
