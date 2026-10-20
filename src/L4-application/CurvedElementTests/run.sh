#!/usr/bin/env bash
# Builds CurvedElementTests and runs it via its own Main (see CurvedElementTestsMain.cs).
#
# Usage: ./run.sh [-c Debug|Release] [-n <MPI processes>] [--no-build] [-- <args for the executable>]
#   -c          build configuration (default: Debug)
#   -n          number of MPI processes; if > 1, mpiexec is used (default: 1)
#   --no-build  skip the build step
set -euo pipefail

CONFIG=Debug
NPROCS=1
BUILD=1
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c) CONFIG="$2"; shift 2 ;;
        -n) NPROCS="$2"; shift 2 ;;
        --no-build) BUILD=0; shift ;;
        -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
        --) shift; EXTRA_ARGS=("$@"); break ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="$SCRIPT_DIR/bin/$CONFIG/net10.0"

if [[ $BUILD -eq 1 ]]; then
    echo "=== building CurvedElementTests ($CONFIG) ==="
    dotnet build "$SCRIPT_DIR/CurvedElementTests.csproj" -c "$CONFIG" -v q
fi

if [[ ! -f "$BIN_DIR/CurvedElementTests.dll" ]]; then
    echo "CurvedElementTests.dll not found in $BIN_DIR; build first." >&2
    exit 1
fi

echo "=== running CurvedElementTests ($CONFIG, $NPROCS process(es)) ==="
cd "$BIN_DIR"
STATUS=0
if [[ $NPROCS -gt 1 ]]; then
    mpiexec -n "$NPROCS" dotnet CurvedElementTests.dll "${EXTRA_ARGS[@]}" || STATUS=$?
else
    dotnet CurvedElementTests.dll "${EXTRA_ARGS[@]}" || STATUS=$?
fi

if [[ $STATUS -eq 0 ]]; then
    echo "=== PASSED ==="
else
    echo "=== FAILED (exit code $STATUS) ===" >&2
fi
exit $STATUS
