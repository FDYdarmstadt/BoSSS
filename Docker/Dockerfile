FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG USERNAME=bosss
ARG UID=1000
ARG GID=1000

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		ca-certificates \
		curl \
		git \
		gnupg \
		libopenmpi-dev \
		openmpi-bin \
		octave \
		gnuplot \
		openssh-client \
		sudo \
		wget \
		nano \
	&& rm -rf /var/lib/apt/lists/*

# Add Microsoft package feed for .NET SDK 8 and 10.
RUN wget https://packages.microsoft.com/config/ubuntu/24.04/packages-microsoft-prod.deb -O /tmp/packages-microsoft-prod.deb \
	&& dpkg -i /tmp/packages-microsoft-prod.deb \
	&& rm /tmp/packages-microsoft-prod.deb \
	&& apt-get update \
	&& apt-get install -y --no-install-recommends \
		dotnet-sdk-8.0 \
		dotnet-sdk-10.0 \
	&& rm -rf /var/lib/apt/lists/*

# Install Node.js (for CLI tools) from NodeSource.
RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash - \
	&& apt-get update \
	&& apt-get install -y --no-install-recommends nodejs \
	&& rm -rf /var/lib/apt/lists/*

# Try to install Claude and Codex CLIs.
# Build does not fail if upstream package names or availability change.
RUN npm install -g @anthropic-ai/claude-code || true
RUN npm install -g @openai/codex || true

# Non-root user by default (recommended for OpenMPI usage and safer dev workflows).
RUN groupadd --gid ${GID} ${USERNAME} \
	&& useradd --uid ${UID} --gid ${GID} -m -s /bin/bash ${USERNAME} \
	&& usermod -aG sudo ${USERNAME} \
	&& echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME} \
	&& chmod 0440 /etc/sudoers.d/${USERNAME}

RUN install -d -o ${USERNAME} -g ${USERNAME} /home/${USERNAME}/.BoSSS/etc \
	&& cat > /home/${USERNAME}/.BoSSS/etc/MatlabConnectorConfig.json <<'EOF'
{
  "Flav": "Octave"
}
EOF

RUN chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/.BoSSS/

RUN install -d -o ${USERNAME} -g ${USERNAME} /home/${USERNAME}/.codex \
	&& cat > /home/${USERNAME}/.codex/config.toml <<'EOF'
approval_policy = "never"
sandbox_mode = "danger-full-access"

model = "gpt-5.6-luna"
model_reasoning_effort = "medium"
service_tier = "default"

[projects."/workspace"]
trust_level = "trusted"

EOF

RUN chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/.codex


RUN install -d -o ${USERNAME} -g ${USERNAME} /home/${USERNAME}/.claude \
	&& cat > /home/${USERNAME}/.claude/settings.json <<'EOF'
{
  "permissions": {
    "defaultMode": "bypassPermissions",
    "deny": [
      "Bash(rm *)",
      "Bash(rmdir *)",
      "Bash(unlink *)",
      "Bash(shred *)",
      "Bash(find * -delete*)",
      "Bash(find * -exec rm*)",
      "Bash(git clean*)",
      "Bash(git rm*)"
    ]
  },
  "skipDangerousModePermissionPrompt": true,
  "model": "opus",
  "effortLevel": "high"
}

EOF

RUN chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/.claude


# Runtime location for full BoSSS host install copied from host-mounted folder.
RUN mkdir -p /home/${USERNAME}/BoSSS-Install \
	&& chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/BoSSS-Install

RUN cat <<'EOF' > /usr/local/bin/bosss-entrypoint.sh
#!/usr/bin/env bash
set -uo pipefail

warn() {
	echo "[bosss-entrypoint] WARNING: $*" >&2
}

configure_git_identity_from_host() {
	local host_gitconfig="${HOST_GITCONFIG:-}"
	local candidate
	local host_name=""
	local host_email=""

	if [[ -z "$host_gitconfig" ]]; then
		for candidate in \
			"/host-home/.gitconfig" \
			"/host/.gitconfig" \
			"/mnt/host-home/.gitconfig" \
			"/run/host-home/.gitconfig"; do
			if [[ -f "$candidate" ]]; then
				host_gitconfig="$candidate"
				break
			fi
		done
	fi

	if [[ -n "$host_gitconfig" && -f "$host_gitconfig" ]]; then
		host_name="$(git config --file "$host_gitconfig" --get user.name || true)"
		host_email="$(git config --file "$host_gitconfig" --get user.email || true)"

		if [[ -n "$host_name" ]]; then
			git config --global user.name "$host_name" || warn "Failed to set git user.name from host config"
		fi

		if [[ -n "$host_email" ]]; then
			git config --global user.email "$host_email" || warn "Failed to set git user.email from host config"
		fi
	fi
}

export BOSSS_INSTALL="${BOSSS_INSTALL:-${HOME}/BoSSS-Install}"
SYNC_MARKER="${HOME}/.BoSSS/.host-sync-done"

configure_git_identity_from_host

# If host BoSSS install is mounted to /host-bosss-install, copy it into user home.
# Copy is done once by default; set BOSSS_FORCE_SYNC=1 to overwrite on every start.
if [[ -d /host-bosss-install ]]; then
	mkdir -p "$BOSSS_INSTALL" || warn "Could not create BOSSS_INSTALL directory: $BOSSS_INSTALL"
	mkdir -p "$(dirname "$SYNC_MARKER")" || warn "Could not create sync-marker directory"
	if [[ "${BOSSS_FORCE_SYNC:-0}" == "1" || ! -f "$SYNC_MARKER" ]]; then
		rm -rf "$BOSSS_INSTALL" || warn "Could not remove existing BOSSS_INSTALL directory"
		mkdir -p "$BOSSS_INSTALL" || warn "Could not recreate BOSSS_INSTALL directory"
		if ! cp -r /host-bosss-install/. "$BOSSS_INSTALL"/; then
			warn "Failed to copy /host-bosss-install into $BOSSS_INSTALL; continuing without sync"
		else
			chown -R "$(id -u)":"$(id -g)" "$BOSSS_INSTALL" || warn "Could not normalize ownership in $BOSSS_INSTALL"
			date -u +"%Y-%m-%dT%H:%M:%SZ" > "$SYNC_MARKER" || warn "Could not update sync marker: $SYNC_MARKER"
		fi
	fi
fi

if [[ -d "$BOSSS_INSTALL/bin/native/linux/amd64-openmpi" ]]; then
	export LD_LIBRARY_PATH="$BOSSS_INSTALL/bin/native/linux/amd64-openmpi:${LD_LIBRARY_PATH:-}"
else
	warn "Native library path not found: $BOSSS_INSTALL/bin/native/linux/amd64-openmpi"
	export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
fi

exec "$@"
EOF

RUN chmod +x /usr/local/bin/bosss-entrypoint.sh
RUN sed -i 's/\r$//' /usr/local/bin/bosss-entrypoint.sh

ENV BOSSS_INSTALL=/home/${USERNAME}/BoSSS-Install
ENV LD_LIBRARY_PATH=/home/${USERNAME}/BoSSS-Install/bin/native/linux/amd64-openmpi

WORKDIR /workspace
USER ${USERNAME}

ENTRYPOINT ["/usr/local/bin/bosss-entrypoint.sh"]
CMD ["bash"]
