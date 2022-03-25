# This Dockerfile builds the React client and API together
FROM continuumio/miniconda3
COPY . /biosnicar
WORKDIR /biosnicar
RUN ls
RUN conda env create -f ./app/api/environment.yml
# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "biosnicar", "/bin/bash", "-c"]


# Build step #1: build the React front end
FROM node:16
WORKDIR /biosnicar
RUN ls
COPY ./app/package.json ./app/

WORKDIR /biosnicar/app
RUN yarn install
RUN yarn build


# Build step #2: build the API with the client as static files
FROM python:3.6
RUN pwd

ENV FLASK_ENV production

EXPOSE 3000
WORKDIR ./app/api
CMD ["gunicorn", "-b", ":3000", "app:app"]