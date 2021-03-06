qsip Package
================
Bram Stone

-   [Background](#background)
    -   [Terminology](#terminology)
-   [Installation](#installation)
-   [Data preparation](#data-preparation)
    -   [Creation of phyloseq object](#creation-of-phyloseq-object)
-   [Creation of phylosip object](#creation-of-phylosip-object)
    -   [Plot enrichment curves](#plot-enrichment-curves)
-   [Specify frequency filtering criteria](#specify-frequency-filtering-criteria)
-   [Calculate enrichment and growth](#calculate-enrichment-and-growth)
-   [Extract data frames from phylosip object](#extract-data-frames-from-phylosip-object)
-   [References](#references)

Background
----------

**Quantitative stable isotope probing (qSIP)** is the combination of stable isotope probing – a foundational technique in the study of ecosystems – with targeted amplicon sequencing data of microbial communities. In conventional stable isotope probing (SIP) experiments, identification of the amount of isotopic enrichment of nucleic acids was done qualitatively – through visual identification and categorization of nucleic acids into either "heavy" or "light" regions. The qSIP approach is to divide a single sample into many different fractions (without *a priori* categorization) along a gradient of increasing densities, and to estimate the shift in density of an individual microbial taxon's nucleic acids based on it's abundance across the many fractions (Hungate *et al.* 2015). Further details on the qSIP methodology may be found in Purcell *et al.* (2019). The core calculation produces estimates of every microbial taxon's proportion of enrichment which are often used as stand-ins for growth. However, **qsip** can also estimate population *per-capita* rates of growth and death (Koch *et al.* 2018).

The **qsip** package is built on the data structures of the **phyloseq** [package](http://joey711.github.io/phyloseq/index.html). **Phyloseq** uses S4 object classification to organize different aspects of microbial sequencing data into a single object (see Data Preparation, below). The purpose of **phyloseq** is to minimize the code necessary to perform common operations on microbial community datasets such as filtering out, or retaining, certain taxonomic lineages or groups of samples.

The **qsip** package supports stable isotope experiments using <sup>18</sup>O, <sup>13</sup>C, and <sup>15</sup>N (Morrissey *et al.* 2018). It is also agnostic towards the amplicon being sequenced. It will support 16S, 18S, ITS, or other amplicon sequencing data.

### Terminology

Conducting a qSIP experiment is, in many ways, an exercise in data organization. In an amplicon sequencing study, a single sequencing sample is usually produced from DNA extracted from a single point of collection (unless samples are pooled, in which case many points of collection yield one sequencing sample). In a qSIP experiment, DNA from each sample is divided into usually more than a dozen fractions which must all be sequenced separately. Because of this, as well as to make this package as easy to use as possible, consistent terminology should be applied to any qSIP experiment.

-   **Replicate** - A single physical sample that has been divided into multiple fractions (usually 10–20). The term *replicate* is used rather than *sample* to avoid confusion during the sequencing preparation process, in which each qSIP fraction must be treated as a separate sample. Furthermore, the term *replicate* makes it clear that these represent the unit of statistical replication and power.
-   **Fraction** - A subsample produced by fractionating a single replicate based on vertical stratification following ultra-centrifugation. *Each fraction must have its own qPCR amplicon abundance measurement as well as a density measurement.*
-   **Feature Table** - Often called an OTU table, ASV table, or species abundance table. A table of abundances or relative frequences for each microbial taxon in each sample. For qSIP analyses, feature tables must include the abundances/frequencies in *each* fraction. For qSIP analyses, this should not be a binary presence-absence table.
-   **Labeled** - A replicate that has been treated with a stable isotope (<sup>18</sup>O, <sup>13</sup>C, or <sup>15</sup>N). Often used to differentiate labeled from unlabeled values in output from the qsip package (as well as represented in equations throughout the published literature).
-   **Light** - A replicate that has **not** been treated with a stable isotope. An unlabeled replicate.
-   **Atom excess fraction (AEF)** - The proportion \[0–1\] of an organism's nucleic acids that are enriched by a stable isotope. Identical to **atom percent excess (APE)** or **atom excess percent (AEP)**.

Installation
------------

**qsip** is not currently on CRAN. The only way to install **qsip** is through Github. See here if you encounter issues [installing phyloseq](http://joey711.github.io/phyloseq/install.html).

``` r
# install devtools and BiocManager
install.packages('devtools')
install.packages('BiocManager')

# install phyloseq and qsip using the utilities on devtools and BiocManagerK
BiocManager::install('phyloseq')
devtools::install_github('bramstone/qsip')

library(phyloseq)
library(qsip)
```

Data preparation
----------------

Here, go into how to organize and prepare experimental data, the feature table, and taxonomic data for phyloseq combination. Spend the most time on experimental data.

### Creation of phyloseq object

Feature tables must be combined with taxonomic and experimental data into a single [phyloseq object](http://joey711.github.io/phyloseq/import-data.html).

Creation of phylosip object
---------------------------

``` r
dat <- specify_qsip(dat,
                    density='Density.g.ml',
                    abund='avg_16S_g_soil',
                    rep_id='sampleID',
                    rep_group='ecosystem',
                    iso='18O',
                    iso_trt='isotope',
                    timepoint='timepoint')

# high-level qSIP-related data will display upon printing
dat
```

### Plot enrichment curves

One of the first diagnostic plots that should be generated from an enrichment experiment is the enrichment curve. Often, this is done even before sequencing efforts begin. The basic idea of the enrichment curve is to plot the change in density of DNA due to the incorporation of stable isotopes.

*With qsip, enrichment curves may be plotted from phyloseq, phylosip, or basic data.frame objects.* Because phylosip objects already have replicates, fractions, densities, and qPCR abundances specified, no other information is necessary. Using phyloseq or data.frame objects requires the specification of this information.

``` r
# plot from phylosip object
plot_curve(dat)

# plot from phyloseq object
plot_curve(phlyo_dat, density='density_column', abund='qPCR_column', iso_trt='isotope_treatment_column')

# plot from data.frame object
plot_curve(exper_dat, density='density_column', abund='qPCR_column', iso_trt='isotope_treatment_column')
```

Specify frequency filtering criteria
------------------------------------

``` r
dat@qsip@filter_levels <- create_filters(2, 5, soft=1)
```

Calculate enrichment and growth
-------------------------------

``` r
# atom excess fraction (AEF) or atom percent excess (APE)
dat <- calc_excess(dat,
                   separate_label=T,
                   filter=T,
                   correction=T)

# per-capita growth
dat <- calc_pop(dat,
                separate_label=T,
                filter=T,
                correction=T,
                growth_model='exponential')

# show new data has been added
dat
names(dat@qsip)
```

Extract data frames from phylosip object
----------------------------------------

``` r
# extract data related to isotopic "excess" and measures from "label"-ed samples
qsip_enrich <- qsmelt(dat, include='excess|label', regex=T)
head(qsip_enrich)

# extract data related to birth and death
qsip_growth <- qsmelt(dat, include='birth|death', abundance=T, regex=T)
head(qsip_growth)
```

References
----------

Hungate BA, Mau RL, Schwartz E *et al.* Quantitative Microbial Ecology through Stable Isotope Probing. *Applied and Environmental Microbiology* 2015;**81**, DOI: [10.1128/AEM.02280-15](https://doi.org/10.1128/AEM.02280-15).

Koch BJ, McHugh TA, Hayer M *et al.* Estimating taxon-specific population dynamics in diverse microbial communities. *Ecosphere* 2018;**9**, DOI: [10.1002/ecs2.2090](https://doi.org/10.1002/ecs2.2090).

Morrissey EM, Mau RL, Schwartz E *et al.* Taxonomic patterns in the nitrogen assimilation of soil prokaryotes. *Environmental Microbiology* 2018;**20**, DOI: [10.1111/1462-2920.14051](https://doi.org/10.1111/1462-2920.14051).

Purcell AM, Dijkstra P, Finley B *et al.* Quantitative Stable Isotope Probing with H218O to Measure Taxon-Specific Microbial Growth. *Methods of Soil Analysis* 2019;**4**, DOI: [10.2136/msa2018.0083](https://doi.org/10.2136/msa2018.0083).
