FROM ghcr.io/tgoelles/python_docker:v0.3.5_py3.9


ENV DEBIAN_FRONTEND=noninteractive


COPY .devcontainer/environment.yml /tmp/conda-tmp/
RUN /opt/conda/bin/conda env update -n base -f /tmp/conda-tmp/environment.yml
