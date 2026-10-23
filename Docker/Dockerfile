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
# Base images may already contain the requested IDs (e.g. ubuntu at 1000:1000).
# Reuse the group and rename/move an existing user rather than duplicating IDs.
RUN test "${UID}" -ne 0 \
	&& if ! getent group "${GID}" > /dev/null; then \
		groupadd --gid "${GID}" "${USERNAME}"; \
	fi \
	&& if getent passwd "${UID}" > /dev/null; then \
		existing_user="$(getent passwd "${UID}" | cut -d: -f1)"; \
		if [ "$existing_user" != "${USERNAME}" ]; then \
			usermod --login "${USERNAME}" "$existing_user" || exit 1; \
		fi; \
		usermod --home "/home/${USERNAME}" --move-home \
			--gid "${GID}" --shell /bin/bash "${USERNAME}"; \
	else \
		useradd --uid "${UID}" --gid "${GID}" -m -s /bin/bash "${USERNAME}"; \
	fi \
	&& usermod -aG sudo "${USERNAME}" \
	&& echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME} \
	&& chmod 0440 /etc/sudoers.d/${USERNAME}

COPY public/Docker/MatlabConnectorConfig.json /home/${USERNAME}/.BoSSS/etc/MatlabConnectorConfig.json

RUN chown -R ${USERNAME}:${GID} /home/${USERNAME}/.BoSSS/

COPY public/Docker/config.toml /home/${USERNAME}/.codex/config.toml

RUN chown -R ${USERNAME}:${GID} /home/${USERNAME}/.codex


COPY public/Docker/settings.json /home/${USERNAME}/.claude/settings.json

RUN chown -R ${USERNAME}:${GID} /home/${USERNAME}/.claude


# Runtime location for full BoSSS host install copied from host-mounted folder.
RUN mkdir -p /home/${USERNAME}/BoSSS-Install \
	&& chown -R ${USERNAME}:${GID} /home/${USERNAME}/BoSSS-Install

# Copy scripts/configuration explicitly: legacy Docker builders can produce empty heredocs.
COPY public/Docker/bosss-entrypoint.sh /usr/local/bin/bosss-entrypoint.sh

RUN chmod +x /usr/local/bin/bosss-entrypoint.sh

ENV BOSSS_INSTALL=/home/${USERNAME}/BoSSS-Install
ENV LD_LIBRARY_PATH=/home/${USERNAME}/BoSSS-Install/bin/native/linux/amd64-openmpi

WORKDIR /workspace
USER ${USERNAME}

ENTRYPOINT ["/usr/local/bin/bosss-entrypoint.sh"]
CMD ["bash"]
