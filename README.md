---
author: JR Ferrer-Paris
date: 28 Oct 2026
---

# Cross-reference of the VSNT to the IUCN GET

Source classification is the *Vegetation survey of the Northern Territory, Australia* (VSNT).

Target classification is the *IUCN Global Ecosystem Typology* (IUCN GET).

Cross-referencing completed by MG Tozer, following IUCN (2025).

Workflow for synthesis and outputs created by JR Ferrer-Paris using the *ecoSynthesis* and *fasterize* packages.

## Required data

The VSNT has two components:

- [Description of map units](https://hdl.handle.net/10070/674164)
- [Geospatial data](https://ftp-dlrm.nt.gov.au/main.html?download&weblink=e78bd012813d86c84c37681c72021704&realfilename=VegResources_ntveg_1M.zip) (Link provided by MGT, metadata is missing)

For the IUCN GET we are using the version provided by the *ecoSynthesis* package, which is based on the [IUCN GET vocabulary](https://w3id.org/iucn-get/). Note that MGT cross-references include a draft unit *T5.6* which is not yet included in the vocabulary.

Cross-referencing table in paired column format were provided by MGT and are available in folder [`xrefdata`](/xrefdata/). The files include header with metadata.

## Software notes

Main code in folder [`scripts`](/scripts/) is written in [R](https://cran.r-project.org/) and uses several packages from the [tidyverse](https://tidyverse.org/packages/), [rspatial](https://rspatial.org/index.html) and [r-spatial](https://r-spatial.org/) collections.

For working with cross-referencing data we use the development version of [*ecoSynthesis*](https://github.com/Ecosystem-Indicators-Workflows/ecoSynthesis). The [*ecoSynthesis* documentation](https://ecosystem-indicators-workflows.github.io/ecoSynthesis/) is work in progress.

Install from GitHub:

```
pak::pak("Ecosystem-Indicators-Workflows/ecoSynthesis")
```

We use [*fasterize*](https://hypertidy.github.io/fasterize/) for fast polygon-to-raster conversion, burn polygon shapes and/or values into pixels. Install the development version of fasterize with 

```
pak::pak("hypertidy/fasterize")
```

*fasterize* uses *Rcpp* and thus requires a compile toolchain to install from source. Testing (and for normal use of *sf* objects) requires *sf*, which requires *GDAL*, *GEOS*, and *PROJ* to be installed.

## Related resources

This repository is part of the **Ecosystem Indicator Workflows** project. Visit the [project page](https://ardc.edu.au/project/ecosystem-indicator-workflows/) (ARDC) and [project website](https://ecosystem-indicators-workflows.github.io/).

## Acknowledgements

The Ecosystem Indicator Workflows is a co-investment partnership with the Australian Research Data Commons (ARDC) through the Planet Research Data Commons (DOI: [10.3565/j5sq-a757](http://doi.org/10.3565/j5sq-a757)). The ARDC is enabled by the Australian Government’s National Collaborative Research Infrastructure Strategy (NCRIS).  

## References

Ferrer-Paris J, Gould E, Singh A (2026). ecoSynthesis: Synthesis of cross-referencing information between two classification systems. R package version 0.0.0.9002, https://github.com/Ecosystem-Indicators-Workflows/ecoSynthesis.

IUCN. 2025. Standards, Methods and Guidelines for Cross-Referencing Ecosystem Classifications and Maps to the IUCN Global Ecosystem Typology (version Version 1.0). Edited by D. A. Keith, A. Dayaram, A. Driver, J. R. Ferrer-Paris, R. Francis, S. Hay, N. J. Murray, et al. Gland, Switzerland: IUCN. https://doi.org/10.2305/CPPC6126.

Wilson, B. A.. (1990, February 26) Vegetation survey of the Northern Territory, Australia 1990 - explanatory notes to accompany 1: 1000000 map sheets. Retrieved 2026, October 27, from https://hdl.handle.net/10070/674164.