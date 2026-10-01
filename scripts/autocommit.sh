#!/usr/bin/env bash
# Automatic commit + push of the project repository (used by the Claude Code Stop hook and by long-run orchestrators).
# Commits only if there are changes (respecting .gitignore); never fails the caller.
# Usage: scripts/autocommit.sh ["message"]
cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/..}" || exit 0
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || exit 0
git add -A >/dev/null 2>&1
if git diff --cached --quiet; then
  exit 0
fi
msg="${1:-Auto-commit: $(date '+%Y-%m-%d %H:%M %Z')}"
git commit -q -m "$msg

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" >/dev/null 2>&1 || exit 0
git push -q >/dev/null 2>&1 || echo "autocommit: push failed (commit kept locally)" >&2
exit 0
