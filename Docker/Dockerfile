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
default_permissions = "danger-full-access"
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
      "Bash(git commit*)",
      "Bash(git * commit*)",
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
set -euo pipefail

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
EOF

RUN chmod +x /usr/local/bin/bosss-entrypoint.sh

ENV BOSSS_INSTALL=/home/${USERNAME}/BoSSS-Install
ENV LD_LIBRARY_PATH=/home/${USERNAME}/BoSSS-Install/bin/native/linux/amd64-openmpi

WORKDIR /workspace
USER ${USERNAME}

ENTRYPOINT ["/usr/local/bin/bosss-entrypoint.sh"]
CMD ["bash"]
