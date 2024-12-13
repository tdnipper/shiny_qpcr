# Shiny_qpcr

A shiny app to calculate ddCT from qPCR results. Delivered via docker/podman container.

## Overview

This app is designed to take qPCR CT data and return a downloadable dataframe of ddCT values and a graph to inspect the data by group. The app is packaged in a container that can be run by docker or podman and contains all dependencies. The dockerfile is included in the root directory of the container.

## Input data

Input should be formatted in a .csv, .xls, or .xlsx file in three columns:
| Sample Name | Target Name | CT |
| ------------------- | ------------------ | ---|
| Group 1 | Gene 1 | 20 |
| Group 1 | Gene 1 | 21 | 
| Group 1 | Gene 1 | 19 |
| Group 1 | Gene 2 | 6 |
| Group 1 | Gene 2 | 4 | 
| Group 1 | Gene 2 | 3 | 

...etc. 

Samples with the same Group and Gene will be averaged together to take the mean of the CT. If you want biological replicates to be handled seperately, append a suffix to each group name to keep them apart (e.g. WT_1, WT_2).

## Running the container

The container is hosted at `ghcr.io/tdnipper/shiny_qpcr` and can be run via [docker](https://www.docker.com/) or [podman](https://podman.io/get-started).
When running the container, expose port 8000 so that the web UI can be accessed on the local host. For example: `podman run --rm -p 8000:8000 ghcr.io/tdnipper/shiny_qpcr`. Running this via CLI will generate a link to `localhost:8000` that accesses the web UI. If you're using wsl, make sure to navigate to `localhost:8000` instead of `0.0.0.0:8000`.

## Running natively

If you want to run this without using a container, you can run the app natively. Clone the repository, create a python virtual environment, and use requirements.txt to install dependencies. Then activate the venv and run the program using `shiny run app.py`. 