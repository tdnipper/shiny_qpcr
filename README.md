# Shiny_qpcr

A shiny app to calculate ddCT from qPCR results. Delivered via docker/podman container.

## Overview

This app is designed to take qPCR CT data and return a downloadable dataframe of ddCT values and a graph to inspect the data by group. The app is packaged in a container that can be run by docker or podman and contains all dependencies. The dockerfile is included in the root directory of the container.

## Input data

Input should be formatted in a .csv, .xls, or .xlsx file as:
| Sample (Group) Name | Target (Gene) Name | CT |
| ------------------- | ------------------ | ---|
| Group 1 | Gene 1 | 20 |

## Running the container

The container is hosted at `quay.io/tdnipper/shiny_qpcr` and can be run via [docker](https://www.docker.com/) or [podman](https://podman.io/get-started).
When running the container, expose port 8000 so that the web UI can be accessed on the local host. For example: `podman run --rm -p 8000:8000 quay.io/tdnipper/shiny_qpcr`. Running this via CLI will generate a link to `localhost:8000` that accesses the web UI.