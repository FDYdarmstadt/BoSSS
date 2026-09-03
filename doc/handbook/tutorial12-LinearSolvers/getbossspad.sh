#!/usr/bin/env bash
set -eu

dotnet publish ../../../src/L4-application/BoSSSpad/BoSSSpad.csproj -c Release -v q -o .
