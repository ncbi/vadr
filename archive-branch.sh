#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <existing-branch> [remote]" >&2
  echo "Example: $0 bug-cdsstopn-seqend origin" >&2
  exit 1
fi

src_branch="$1"
remote="${2:-origin}"
archive_name="$src_branch"
archive_name="${archive_name#ABD-}"
dst_branch="z-archive/${archive_name}"

if [[ "$src_branch" == "$dst_branch" ]]; then
  echo "Error: source branch already has z-archive/ prefix: $src_branch" >&2
  exit 1
fi

if ! git rev-parse --git-dir >/dev/null 2>&1; then
  echo "Error: not inside a git repository" >&2
  exit 1
fi

if ! git show-ref --verify --quiet "refs/heads/$src_branch"; then
  echo "Error: local branch does not exist: $src_branch" >&2
  exit 1
fi

if git show-ref --verify --quiet "refs/heads/$dst_branch"; then
  echo "Error: destination local branch already exists: $dst_branch" >&2
  exit 1
fi

if git ls-remote --exit-code --heads "$remote" "$dst_branch" >/dev/null 2>&1; then
  echo "Error: destination remote branch already exists: $remote/$dst_branch" >&2
  exit 1
fi

echo "Step 1/4: rename local branch: $src_branch -> $dst_branch"
git branch -m "$src_branch" "$dst_branch"

echo "Step 2/4: push new archived branch to $remote"
git push "$remote" "$dst_branch"

echo "Step 3/4: delete old branch on $remote (if it exists)"
if git ls-remote --exit-code --heads "$remote" "$src_branch" >/dev/null 2>&1; then
  git push "$remote" --delete "$src_branch"
else
  echo "  Remote branch $remote/$src_branch does not exist; skipping delete"
fi

echo "Step 4/4: set upstream for local archived branch"
git push -u "$remote" "$dst_branch"

echo "Done: $src_branch archived as $dst_branch"
