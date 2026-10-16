#!/usr/bin/env python3
"""Flag added lines that re-roll something pic already provides.

Checks only the lines a change ADDS (git diff against a base), in src/ and
app/, so existing code does not fail CI; the pr-review skill (section 10)
covers the judgement calls this cannot make.

    python3 tools/pic_reuse_check.py origin/main

Flags, with the pic replacement (docs/pic_api.md lists all of pic):
  - system_clock / cpu_time / date_and_time / omp_get_wtime
        -> pic_timer (timer_type) or pic_clock (monotonic_ms/us)
  - a new integer-to-string function (itoa, int2str, to_str*, ...)
        -> pic_strings (to_string)
  - a new sort routine (sort*, *_sort)
        -> pic_sorting (sort, sort_index, ord_sort, radix_sort)
  - a new mean / variance / std_dev function
        -> pic_statistics

Device code is exempt (pic is host code): end the flagged line with
`! pic-reuse: device` or `! pic-reuse: <reason>` and say why in the PR.
"""
import re
import subprocess
import sys

RULES = [
    (re.compile(r"\b(system_clock|cpu_time|date_and_time|omp_get_wtime)\b", re.I),
     "time with pic_timer (timer_type) or pic_clock (monotonic_ms/us)"),
    (re.compile(r"\bfunction\s+(itoa|i2s|int2str|num2str|str_of\w*|to_str\w*)\b", re.I),
     "convert numbers to text with pic_strings' to_string"),
    (re.compile(r"\b(subroutine|function)\s+(sort\w*|\w+_sort)\b", re.I),
     "sort with pic_sorting (sort, sort_index, ord_sort, radix_sort)"),
    (re.compile(r"\bfunction\s+(mean|variance|std_?dev)\b", re.I),
     "use pic_statistics (mean, variance, std_dev, min_max)"),
]
WAIVER = re.compile(r"!\s*pic-reuse:\s*\S", re.I)
PATHS = re.compile(r"^(src|app)/.*\.[fF]90$")


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else "origin/main"
    diff = subprocess.run(["git", "diff", "-U0", f"{base}...HEAD", "--", "src", "app"],
                          capture_output=True, text=True, check=True).stdout
    path, line_no, hits = None, 0, 0
    for line in diff.splitlines():
        if line.startswith("+++ "):
            p = line[6:] if line.startswith("+++ b/") else None
            path = p if p and PATHS.match(p) and "/generated/" not in p else None
            continue
        m = re.match(r"@@ -\S+ \+(\d+)", line)
        if m:
            line_no = int(m.group(1))
            continue
        if not line.startswith("+") or line.startswith("+++"):
            continue
        text = line[1:]
        if path:
            code = text.split("!", 1)[0]
            for rx, advice in RULES:
                if rx.search(code) and not WAIVER.search(text):
                    print(f"::error file={path},line={line_no}::pic provides this: {advice} "
                          "(or waive with `! pic-reuse: <reason>`, e.g. device code)")
                    hits += 1
        line_no += 1
    if hits:
        print(f"pic_reuse_check: {hits} added line(s) re-roll something pic provides; see docs/pic_api.md")
        return 1
    print("pic_reuse_check: no added line re-rolls a pic facility")
    return 0


if __name__ == "__main__":
    sys.exit(main())
