# BoSSS Dev Docker Image

This folder provides a Docker image for BoSSS development on Ubuntu 24.04.

Included tools:
- git
- OpenMPI
- .NET SDK 8.0
- .NET SDK 10.0
- Node.js
- best-effort installation of Claude and Codex CLIs

## Build

Run from repository root:
```bash
docker build -t bosss-dev -f public/Docker/Dockerfile \
  --build-arg UID=$(id -u) --build-arg GID=$(id -g) .
```

## Run with current repo mounted
```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -w /workspace \
  bosss-dev
```
## Run with host BOSSS_INSTALL mounted and copied into container user home

Assuming the host variable BOSSS_INSTALL points to the full BoSSS install tree:
```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -v "$BOSSS_INSTALL":/host-bosss-install:ro \
  -w /workspace \
  bosss-dev
```

At startup, the entrypoint copies /host-bosss-install to $HOME/BoSSS-Install once, then exports BOSSS_INSTALL and updates LD_LIBRARY_PATH.

To force re-sync on every start:
```bash
docker run --rm -it \
  -e BOSSS_FORCE_SYNC=1 \
  -v "$(pwd)":/workspace \
  -v "$BOSSS_INSTALL":/host-bosss-install:ro \
  -w /workspace \
  bosss-dev
```
## User and MPI note

The container default user is non-root (bosss), which is the recommended setup for OpenMPI.
