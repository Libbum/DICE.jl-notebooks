# To build and run this Dockerfile:
# docker build -t myimage .
# docker run -p 8888:8888 myimage

# This Docker has been tested on x86 and ARM.
# Originally adapted from: https://github.com/jupyter/docker-stacks/blob/master/datascience-notebook/Dockerfile

# Depends on:
# https://github.com/jupyter/docker-stacks/blob/master/base-notebook/Dockerfile
# https://github.com/jupyter/docker-stacks/blob/master/minimal-notebook/Dockerfile

FROM "jupyter/minimal-notebook"

USER root

# Julia installation
# Default values can be overridden at build time
# (ARGS are in lower case to distinguish them from ENV)
# Check https://julialang.org/downloads/
ARG julia_version="1.7.2"

# Julia dependencies
# install Julia packages in /opt/julia instead of ${HOME}
ENV JULIA_DEPOT_PATH=/opt/julia \
    JULIA_PKGDIR=/opt/julia \
    JULIA_VERSION="${julia_version}"

WORKDIR /tmp

# This will download the right architecture for the current system.
RUN set -x && \
    julia_arch=$(uname -m) && \
    julia_short_arch="${julia_arch}" && \
    if [ "${julia_short_arch}" == "x86_64" ]; then \
      julia_short_arch="x64"; \
    fi; \
    julia_installer="julia-${JULIA_VERSION}-linux-${julia_arch}.tar.gz" && \
    julia_major_minor=$(echo "${JULIA_VERSION}" | cut -d. -f 1,2) && \
    mkdir "/opt/julia-${JULIA_VERSION}" && \
    wget -q "https://julialang-s3.julialang.org/bin/linux/${julia_short_arch}/${julia_major_minor}/${julia_installer}" && \
    tar xzf "${julia_installer}" -C "/opt/julia-${JULIA_VERSION}" --strip-components=1 && \
    rm "${julia_installer}" && \
    ln -fs /opt/julia-*/bin/julia /usr/local/bin/julia

# Create JULIA_PKGDIR \
RUN mkdir "${JULIA_PKGDIR}" && \
    chown "${NB_USER}" "${JULIA_PKGDIR}" && \
    fix-permissions "${JULIA_PKGDIR}"

USER $NB_UID

# Add packages and precompile
RUN julia -e 'import Pkg; Pkg.update()'
RUN julia -e 'import Pkg; Pkg.add("IJulia"); using IJulia'

# For some reason, this takes forever.
#RUN fix-permissions /home/$NB_USER

WORKDIR "${HOME}"

COPY Manifest.toml Project.toml ./

RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile();'

RUN rmdir work
COPY 2013 2013
COPY 2016 2016
COPY GAMS GAMS

ENV DOCKER_STACKS_JUPYTER_CMD=lab
