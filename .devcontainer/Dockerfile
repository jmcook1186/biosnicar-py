FROM ghcr.io/tgoelles/python_docker:v0.3.5_py3.8


ENV DEBIAN_FRONTEND=noninteractive


RUN curl -fsSL https://deb.nodesource.com/setup_current.x | bash - \
    && apt-get update \
    # cmdtest incompatible with yarn
    && apt remove cmdtest \
    # add packages here like
    && apt-get -y install nodejs \
    # Clean up
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

RUN /usr/bin/npm install -g yarn


COPY .devcontainer/environment.yml /tmp/conda-tmp/
RUN /opt/conda/bin/conda env update -n base -f /tmp/conda-tmp/environment.yml