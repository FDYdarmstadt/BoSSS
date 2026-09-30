# BoSSS Dev Docker Image

This folder provides a Docker image for BoSSS development on Ubuntu 24.04 (main system) and Ubuntu 26.04 (currently experimental).

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

Windows CMD:
```bat
docker build -t bosss-dev -f public/Docker/Dockerfile .
```

Windows PowerShell:
```powershell
docker build -t bosss-dev -f public/Docker/Dockerfile .
```

Note:
- On Windows, Linux-style `id -u` / `id -g` are usually not available in CMD/PowerShell.
- The image defaults to UID/GID `1000`, which is typically fine for Docker Desktop bind mounts.

## Run with host BOSSS_INSTALL mounted and copied into container user home

Assuming the host variable BOSSS_INSTALL points to the full BoSSS install tree:
```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -v "$BOSSS_INSTALL":/host-bosss-install:ro \
  -w /workspace \
  bosss-dev
```

Windows CMD:
```bat
docker run --rm -it -v "%cd%":/workspace -v "%BOSSS_INSTALL%":/host-bosss-install:ro -w /workspace bosss-dev
```

Windows PowerShell:
```powershell
docker run --rm -it -v "${PWD}:/workspace" -v "${env:BOSSS_INSTALL}:/host-bosss-install:ro" -w /workspace bosss-dev
```

If you do not want to mount a host BoSSS install, omit the second `-v` argument:

Windows CMD:
```bat
docker run --rm -it -v "%cd%":/workspace -w /workspace bosss-dev
```

Windows PowerShell:
```powershell
docker run --rm -it -v "${PWD}:/workspace" -w /workspace bosss-dev
```

At startup, the entrypoint-script (`/usr/local/bin/bosss-entrypoint.sh`) 
copies /host-bosss-install to $HOME/BoSSS-Install once, then exports `BOSSS_INSTALL` and updates `LD_LIBRARY_PATH`.

Implementation detail:
- The one-time sync marker is stored in `$HOME/.BoSSS/.host-sync-done`.
- This avoids permission issues on some Windows-hosted mounts where writing metadata inside the copied install tree can fail.

To force re-sync on every start:
```bash
docker run --rm -it \
  -e BOSSS_FORCE_SYNC=1 \
  -v "$(pwd)":/workspace \
  -v "$BOSSS_INSTALL":/host-bosss-install:ro \
  -w /workspace \
  bosss-dev
```

Windows CMD:
```bat
docker run --rm -it -e BOSSS_FORCE_SYNC=1 -v "%cd%":/workspace -v "%BOSSS_INSTALL%":/host-bosss-install:ro -w /workspace bosss-dev
```

Windows PowerShell:
```powershell
docker run --rm -it -e BOSSS_FORCE_SYNC=1 -v "${PWD}:/workspace" -v "${env:BOSSS_INSTALL}:/host-bosss-install:ro" -w /workspace bosss-dev
```

## Further useful docker commands
* List all images: `docker image ls`
* List all containers: `docker ps -a`
* Delete image: `docker image rm <image_id>`
* Delete container: `docker rm <container_id>`
* Start and re-attach to a stopped container: `docker start -ai <container_id>`

## User and MPI note

The container default user is non-root (bosss), which is the recommended setup for OpenMPI.
