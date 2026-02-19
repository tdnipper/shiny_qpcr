# Shiny_qpcr

A shiny app to analyze qPCR CT data (ddCT, dCT, and enrichment). Delivered via docker/podman container.

## Overview

This app is designed to take qPCR CT data and return downloadable results and plots. Supported analysis modes are ddCT relative expression, dCT fold change between conditions, and RIP/ChIP-qPCR enrichment. The app is packaged in a container that can be run by docker or podman and contains all dependencies. The dockerfile is included in the root directory of the container.

## Input data

The app accepts `.csv`, `.xls`, or `.xlsx` files in one of two formats.

### Flat format (recommended for manual data entry)

Four required columns, plus optional extras:

| Group | Condition | Target Name | CT |
|-------|-----------|-------------|-----|
| WT | untreated | Gene1 | 20.1 |
| WT | untreated | Gene1 | 20.4 |
| WT | treated | Gene1 | 18.3 |
| KO | untreated | Gene1 | 22.0 |

- **Group** — the biological group or genotype (e.g. `WT`, `KO`, `anti_Flag`).
- **Condition** — the experimental condition applied to that group (e.g. `untreated`, `treated`, `input`).
- **Target Name** — the gene or primer target name.
- **CT** — the raw CT value. `Undetermined` and `NTC` are treated as missing and dropped automatically.

An optional **Biological Replicate** column can be included to identify technical replicates. Rows sharing the same Group, Condition, Target Name, and Biological Replicate value are treated as technical replicates and averaged before analysis. This allows, for example, duplicate wells on a plate to be collapsed into a single biological replicate mean.

| Group | Condition | Target Name | CT | Biological Replicate |
|-------|-----------|-------------|-----|----------------------|
| WT | untreated | Gene1 | 20.1 | 1 |
| WT | untreated | Gene1 | 20.4 | 1 |
| WT | untreated | Gene1 | 18.9 | 2 |
| WT | treated | Gene1 | 18.3 | 1 |

In this example, the first two rows are technical replicates of biological replicate 1 and will be averaged to a single CT before analysis. The third row is a distinct biological replicate (replicate 2) and will be kept separate.

**Plots display each biological replicate as an individual dot.** Statistical tests (t-test, ANOVA, etc.) are computed on the biological replicate values — not on the mean across replicates.

### QuantStudio 7 raw export (auto-detected)

Raw `.xlsx` files exported directly from a QuantStudio 7 instrument are accepted without any preprocessing. The app detects these files automatically by the presence of a `Results` sheet.

The **Sample Name** column in the instrument file encodes both Group and Condition, separated by an underscore. The split is always made on the **last** underscore, so group names that themselves contain underscores are handled correctly:

| Sample Name in instrument file | Group parsed | Condition parsed |
|-------------------------------|-------------|-----------------|
| `WT_untreated` | `WT` | `untreated` |
| `KO_treated` | `KO` | `treated` |
| `anti_Flag_enrich` | `anti_Flag` | `enrich` |
| `anti_Flag_input` | `anti_Flag` | `input` |

When naming samples on the instrument, use `GroupName_ConditionName` and ensure every sample name contains at least one underscore. The part after the final underscore becomes the Condition; everything before it becomes the Group.

## Running the container

The container is hosted at `quay.io/tdnipper/shiny_qpcr` and can be run via [docker](https://www.docker.com/) or [podman](https://podman.io/get-started). To run the container, use docker compose or podman-compose to start using the included docker-compose.yml file. This will start the container at http://127.0.0.1:8000.

## Running natively

If you want to run this without using a container, you can run the app natively. Clone the repository, create a python virtual environment, and use requirements.txt to install dependencies. Then activate the venv and run the program using `shiny run app.py`. 
