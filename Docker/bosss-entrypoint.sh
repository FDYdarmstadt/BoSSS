#!/usr/bin/env bash
set -euo pipefail

# Apply only identity values supplied by the host launcher.
if [[ -n "${BOSSS_GIT_USER_NAME:-}" ]]; then
	git config --global user.name "$BOSSS_GIT_USER_NAME"
fi
if [[ -n "${BOSSS_GIT_USER_EMAIL:-}" ]]; then
	git config --global user.email "$BOSSS_GIT_USER_EMAIL"
fi

export BOSSS_INSTALL="${BOSSS_INSTALL:-${HOME}/BoSSS-Install}"

# If host BoSSS install is mounted to /host-bosss-install, copy it into user home.
# Copy is done once by default; set BOSSS_FORCE_SYNC=1 to overwrite on every start.
if [[ -d /host-bosss-install ]]; then
	mkdir -p "$BOSSS_INSTALL"
	if [[ "${BOSSS_FORCE_SYNC:-0}" == "1" || ! -f "$BOSSS_INSTALL/.host-sync-done" ]]; then
		rm -rf "$BOSSS_INSTALL"
		mkdir -p "$BOSSS_INSTALL"
		cp -a /host-bosss-install/. "$BOSSS_INSTALL"/
		date -u +"%Y-%m-%dT%H:%M:%SZ" > "$BOSSS_INSTALL/.host-sync-done"
	fi
fi

export LD_LIBRARY_PATH="$BOSSS_INSTALL/bin/native/linux/amd64-openmpi:${LD_LIBRARY_PATH:-}"

exec "$@"
