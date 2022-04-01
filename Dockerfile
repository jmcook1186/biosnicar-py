# This Dockerfile builds the React client and API together
FROM continuumio/miniconda3 as build-step0
COPY . .

RUN conda env create -f ./app/api/environment.yml
# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "biosnicar", "/bin/bash", "-c"]


# Build step #1: build the React front end
FROM node:16 as build-step1

COPY ./app/package*.json ./app/
WORKDIR ./app
RUN npm install

COPY . .
WORKDIR ./app
RUN npm run build


# Build step #2: build the API with the client as static files
FROM python:3.6 as build-step2
RUN pwd
WORKDIR ./app/api

ENV FLASK_ENV production

EXPOSE 3000
WORKDIR ./app/api
CMD ["gunicorn", "-b", ":3000", "app:app"]