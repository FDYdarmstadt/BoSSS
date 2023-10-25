#!/bin/bash

# Check if git is available
if ! command -v git &> /dev/null; then
    echo "" > commit_hash.txt
    exit 0
fi

#
# Note: all the following is only to make sure that the `commit_hash.txt` is only updated if the commit hash is different;
# otherwise, the incremental build will not work
#

CURRENT_HASH=""
if [ -f "commit_hash.txt" ]; then
    CURRENT_HASH=$(cat commit_hash.txt)
fi

NEW_HASH=$(git rev-parse HEAD)

if [ "$CURRENT_HASH" != "$NEW_HASH" ]; then
    echo $NEW_HASH > commit_hash.txt
fi