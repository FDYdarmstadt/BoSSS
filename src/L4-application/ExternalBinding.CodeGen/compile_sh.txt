#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
# DOTNET_HOST_NATIVE can override SDK host-pack discovery.
dotnet_root=${DOTNET_ROOT:-/usr/share/dotnet}
if [ -z "${DOTNET_HOST_NATIVE:-}" ]; then
    DOTNET_HOST_NATIVE=$(find "$dotnet_root/packs/Microsoft.NETCore.App.Host.linux-x64" -path '*/10.*/runtimes/linux-x64/native' -type d | sort -V | tail -1)
fi
: "${DOTNET_HOST_NATIVE:?Install the .NET 10 SDK native hosting pack}"
"${CXX:-g++}" -std=c++11 -m64 -fPIC -Wall -Wextra -DNETHOST_USE_AS_STATIC \
    -I"$DOTNET_HOST_NATIVE" ./*.cpp "$DOTNET_HOST_NATIVE/libnethost.a" \
    -ldl -pthread -o ExtBindingTest.out
