FROM ghcr.io/tgoelles/python_docker:v0.3.2_py3.8


ENV DEBIAN_FRONTEND=noninteractive

COPY .devcontainer/environment.yml /tmp/conda-tmp/
RUN /opt/conda/bin/conda env update -n base -f /tmp/conda-tmp/environment.yml