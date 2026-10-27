"""Join terco's per-class times (nsys sqlite) with GPU4PySCF's (g4p_jk.py log).

    python join_profile.py terco.sqlite g4p.log [out.json]

terco: the LAST "eri: bins (hgp/rys)" NVTX range (kernels grouped by the
four-digit class in their names) plus the last "mdg: (xx|xx)" range of
each MD exchange class; MD J is reported separately. Classes are put in
angular-momentum-canonical order so both codes' keys match.
"""
import sqlite3, collections, re, sys, json
L = 'spdfg'
def canon(a, b, c, d):
    p = tuple(sorted((a, b), reverse=True)); q = tuple(sorted((c, d), reverse=True))
    if q > p: p, q = q, p
    return '(%s%s|%s%s)' % (L[p[0]], L[p[1]], L[q[0]], L[q[1]])
sq, glog = sys.argv[1], sys.argv[2]
g = collections.Counter()
for l in open(glog):
    m = re.match(r'\((\w)(\w)\|(\w)(\w)\) wall time ([\d.]+)', l)
    if m: g[canon(*[L.index(x) for x in m.groups()[:4]])] += float(m[5])/1e6
c = sqlite3.connect(sq)
last = {}
for n, a, b in c.execute("select coalesce(e.text,s.value), e.start, e.end from NVTX_EVENTS e "
                         "left join StringIds s on s.id=e.textId where e.end is not null order by e.start"):
    last[n] = (a, b)
t = collections.Counter()
a, b = last['eri: bins (hgp/rys)']
for kn, ks, ke in c.execute("select s.value,k.start,k.end from CUPTI_ACTIVITY_KIND_KERNEL k "
                            "join StringIds s on s.id=k.shortName"):
    if ks >= a and ke <= b:
        m = re.search(r'(\d)(\d)(\d)(\d)', kn)
        if m: t[canon(*map(int, m.groups()))] += (ke - ks)/1e9
for n, (s0, s1) in last.items():
    m = re.match(r'mdg: \((\w)(\w)\|(\w)(\w)\)', n)
    if m: t[canon(*[L.index(x) for x in m.groups()])] += (s1 - s0)/1e9
mdj = (last['eri: md J'][1] - last['eri: md J'][0])/1e9 if 'eri: md J' in last else 0.0
agg = collections.defaultdict(lambda: [0.0, 0.0])
for k in set(t) | set(g):
    lm = max(L.index(x) for x in k if x in L)
    agg[L[lm]][0] += t.get(k, 0); agg[L[lm]][1] += g.get(k, 0)
print('max l   terco K   g4p J+K')
for l in L:
    if l in agg: print('  %s   %7.2f   %7.2f' % (l, *agg[l]))
print('terco classes %.2f + MD J %.2f ; g4p %.2f' % (sum(t.values()), mdj, sum(g.values())))
print('largest gaps (terco - g4p):')
for k in sorted(set(t) | set(g), key=lambda k: -(t.get(k, 0) - g.get(k, 0)))[:25]:
    print('  %s terco %.3f g4p %.3f' % (k, t.get(k, 0), g.get(k, 0)))
if len(sys.argv) > 3: json.dump({'terco': t, 'g4p': g, 'mdj': mdj}, open(sys.argv[3], 'w'), indent=0)
