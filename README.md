# Data and code for the manuscript
## Linking Functional Trait-Space Structure to Stream Community Invasion Opportunity across Sub-basins

This repository contains R scripts and data for community-data preparation,
invasibility calculations, and functional and taxonomic diversity metrics.

## Requirements

R 4.5.2 was used by the author. RStudio is optional. No specialised hardware
is required. Install the packages used by Parts 1-3:

```r
install.packages(c("tidyverse", "readxl", "stringr", "sf", "geosphere",
                   "purrr", "mFD", "funrar", "betapart"))
```

## Setup

1. Download the repository and extract `Supplementary files.zip` into the
   project folder.
2. Set the R working directory to that folder.
3. Check the input and output paths in each script. Part 1 reads the annual
   AREMP Excel files named at its beginning and creates `out_watershed/`.
   Part 2 reads the NAS occurrence table, site and watershed attribute
   tables, HUC8 shapefile and watershed-richness table named in the code.
   Keep all components of the shapefile together.
4. Part 2 expects `ALL_YEARS_watershed_richness.csv` in the working directory.
   Copy it from `out_watershed/`, or update that input path.
5. Replace the two personal Windows output paths at the end of Part 2 with
   local paths. Save its invasibility CSV files in the project folder so
   Part 3 can read `Invasibility_WithSpillover.csv`.

A complete rerun requires every input file named in the scripts. Existing
processed tables can be reused when starting at a later stage.

## Run order

Run the scripts from top to bottom in this order:

```r
source("creek analysis part 1 _ exporting files.R")
source("creek analysis part 2 _ invasibility stuff.R")
source("creek analysis part 3_ fd_indices.R")
```

Part 1 exports cleaned community tables and matrices. Part 2 calculates
invasibility scores with and without neighbouring-watershed spillover.
Part 3 calculates diversity indices and joins them to the invasibility data.
Run the accompanying figure scripts afterwards, following their input settings.

## Demonstration and expected outputs

For a shorter demonstration, run Part 3 using these inputs:

- `traits_all.csv`
- `out_watershed/ALL_YEARS_watershed_taxon_matrix_pa_POOLED.csv`
- `Invasibility_WithSpillover.csv`

```r
system.time(source("creek analysis part 3_ fd_indices.R"))
```

Expected outputs are `FD_eigenvalues.csv` and
`FD_indices_with_invasibility.csv`. The latter contains watershed-level
functional-diversity, trait-vacancy, taxonomic beta-diversity and invasibility
metrics. Some functional metrics remain missing for assemblages below their
minimum taxon-richness thresholds.

Installation and execution times depend on the computer and input size.
Numeric benchmarks are not specified here; the command above reports the
demonstration's elapsed time in seconds.

## Using other data

Use the same input columns, identifiers and trait categories as the supplied
tables. Taxon names must match the trait-table keys. Update input paths and
geographic inputs before applying the workflow to another dataset.

## Data sources

Trait-data references are listed in `sources.docx`. Raw sequencing data are
archived under [NCBI BioProject PRJNA1428872](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1428872).

