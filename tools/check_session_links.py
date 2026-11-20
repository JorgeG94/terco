#!/usr/bin/env python3
"""Keep AI-session links out of terco's commits and pull request bodies.

    python3 tools/check_session_links.py BASE HEAD     # commits in BASE..HEAD
    PR_BODY="..." python3 tools/check_session_links.py BASE HEAD

AGENTS.md ("No AI-session links in commits or PR bodies"): a commit message
or PR body must not carry a `Claude-Session:` trailer or a
claude.ai/code/session URL. A `Co-Authored-By:` trailer is fine. Agent tools
add these lines on their own (a commit trailer, a footer on a PR created
through an API), so this check exists to catch the one that slipped through.

Only the commits BASE..HEAD are read, so history already on the base branch
is never re-judged. The PR body comes from the PR_BODY environment variable
(not the command line: a body is untrusted text). Exit status 1 lists each hit.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys

PATTERN = re.compile(r'Claude-Session:|claude\.ai/code/session', re.IGNORECASE)


def commit_hits(base: str, head: str) -> list[str]:
    """`<short sha>: <line>` for each offending line of each commit message in base..head."""
    log = subprocess.run(
        ['git', 'log', '--format=%x00%h%n%B', f'{base}..{head}'],
        check=True, capture_output=True, text=True).stdout
    hits = []
    for chunk in log.split('\0')[1:]:
        sha, _, message = chunk.partition('\n')
        for line in message.splitlines():
            if PATTERN.search(line):
                hits.append(f'commit {sha}: {line.strip()}')
    return hits


def body_hits(body: str) -> list[str]:
    return [f'PR body: {line.strip()}' for line in body.splitlines() if PATTERN.search(line)]


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    hits = commit_hits(argv[0], argv[1]) + body_hits(os.environ.get('PR_BODY', ''))
    for h in hits:
        print(h)
    if hits:
        print(f'check_session_links: {len(hits)} AI-session link(s); remove them '
              '(AGENTS.md: no AI-session links in commits or PR bodies)')
        return 1
    print('check_session_links: no AI-session links')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
