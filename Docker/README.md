# BoSSS Dev Docker Image

This folder provides a Docker image for BoSSS development on Ubuntu 24.04 (main system) and Ubuntu 26.04 (currently experimental).

Included tools:
- git
- OpenMPI
- .NET SDK 8.0
- .NET SDK 10.0
- Node.js
- best-effort installation of Claude and Codex CLIs

## Build and launch

The executable scripts locate the repository relative to their own path, so they
can be called from any working directory.

| Ubuntu | Build only | Launch only | Build then launch | Image |
| --- | --- | --- | --- | --- |
| 24.04 | `./build-ubuntu24-container.sh` | `./launch-ubuntu24-container.sh` | `./buildNlaunch-ubuntu24-container.sh` | `bosss-dev24` |
| 26.04 | `./build-ubuntu26-container.sh` | `./launch-ubuntu26-container.sh` | `./buildNlaunch-ubuntu26-container.sh` | `bosss-dev26` |

Set `BOSSS_INSTALL` before launching:

```bash
export BOSSS_INSTALL=/path/to/BoSSS-Install
./buildNlaunch-ubuntu24-container.sh
```

Each combined script calls the matching build script and then the matching launch
script, stopping on build failure. Launching separately uses the same workspace,
host installation mount and Git identity forwarding.

Both Dockerfiles copy the entrypoint and configuration from files in this directory
so they also work with the legacy Docker builder. Ubuntu 26.04 retains its Codex
settings in `config-ubuntu26.toml`.

## Run with host BOSSS_INSTALL mounted and copied into container user home

Assuming the host variable BOSSS_INSTALL points to the full BoSSS install tree:
```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -v "$BOSSS_INSTALL":/host-bosss-install:ro \
  -w /workspace \
  bosss-dev24
```

At startup, the entrypoint-script (`/usr/local/bin/bosss-entrypoint.sh`) 
copies /host-bosss-install to $HOME/BoSSS-Install once, then exports `BOSSS_INSTALL` and updates `LD_LIBRARY_PATH`.

To force re-sync on every start:
```bash
docker run --rm -it \
  -e BOSSS_FORCE_SYNC=1 \
  -v "$(pwd)":/workspace \
  -v "$BOSSS_INSTALL":/host-bosss-install:ro \
  -w /workspace \
  bosss-dev24
```

## Git identity

The launch scripts read the host's effective `user.name` and `user.email` for this
repository and pass them as `BOSSS_GIT_USER_NAME` and `BOSSS_GIT_USER_EMAIL`.
The shared entrypoint sets those values in the container user's global Git config
on startup. Missing host values leave the corresponding container setting alone.
Only these two identity values are forwarded. Rebuild existing images to include
the updated entrypoint before using the launch-only script.

## Further useful docker commands
* List all images: `docker image ls`
* List all containers: `docker ps -a`
* Delete image: `docker image rm <image_id>`
* Delete container: `docker rm <container_id>`
* Start and re-attach to a stopped container: `docker start -ai <container_id>`

## User and MPI note

The container default user is non-root (bosss), which is the recommended setup for OpenMPI.
