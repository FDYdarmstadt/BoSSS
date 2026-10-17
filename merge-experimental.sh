#!/bin/bash
#
# Forward synchronization: merges the 'public' subtree of the full
# (experimental) repository into the branches 'master' and 'premaster' of this
# (public) repository. Must be called from the public repo.
#
# usage: merge-experimental.sh [REMOTE]
#   REMOTE: URL or path of the experimental repo (optional)
#
# The split commits are created by `git subtree split`, which is deterministic:
# the same history always yields the same hashes. The public repo therefore
# receives the very commits that the experimental repo's own split contains,
# and the backward synchronization (merge-public.sh) does not create duplicates.
#
#   git@git.rwth-aachen.de:bosss1/experimental.git
#   ../BoSSS-experimental   (local test setup)
#
set -e

ROOT_REMOTE=${1:-git@gitlab.hefdi.de:tuda/bosss/bosss-experimental-3.git}
BRANCHES="master premaster"

if ! git diff-index --quiet HEAD --; then
    echo "ERROR: uncommitted changes to tracked files; commit or stash them first." >&2
    exit 1
fi
START_BRANCH=$(git branch --show-current)
echo "Syncing from $ROOT_REMOTE"

for B in $BRANCHES; do
    echo "== Branch $B =="
    ROOT_REF=refs/sync/root/$B
    SPLIT=root-public-$B
    if ! git fetch "$ROOT_REMOTE" "+refs/heads/$B:$ROOT_REF"; then
        echo "   '$B' not found in $ROOT_REMOTE; skipping."
        continue
    fi
    git branch -D "$SPLIT" 2>/dev/null || true
    # git subtree split needs the 'public' prefix in the checked-out tree
    git checkout -q --detach "$ROOT_REF"
    git subtree split -P public -b "$SPLIT" "$ROOT_REF"

    if git show-ref --verify --quiet "refs/heads/$B"; then
        git checkout "$B"
        git merge -X theirs --no-edit -m "Merge branch 'root-public'" "$SPLIT"
    elif git show-ref --verify --quiet "refs/remotes/origin/$B"; then
        git checkout -b "$B" "origin/$B"
        git merge -X theirs --no-edit -m "Merge branch 'root-public'" "$SPLIT"
    else
        echo "   '$B' does not exist in the public repo yet; creating it from the split."
        git branch "$B" "$SPLIT"
    fi
    git branch -D "$SPLIT"
    git update-ref -d "$ROOT_REF"
done

if [ -n "$START_BRANCH" ]; then git checkout "$START_BRANCH"; fi
git branch
