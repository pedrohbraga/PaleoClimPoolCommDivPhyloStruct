# Historical and contemporary processes drive global phylogenetic structure across geographical scales: insights from bat communities

[![badge](https://img.shields.io/static/v1?style=for-the-badge&label=PUBLICATION&message=Open&color=BF616A)](https://doi.org/10.1111/geb.13650) [![badge](https://img.shields.io/static/v1?style=for-the-badge&label=DATA%20at%20DRYAD&message=01&color=B48EAD)](https://doi.org/10.5061/dryad.rjdfn2zgj) [![badge](https://img.shields.io/static/v1?style=for-the-badge&label=CODE%20at%20OSF&message=01&color=8FBCBB)](https://osf.io/amvp5/)

Authors: Pedro Henrique Pereira Braga, Steven Kembel, Pedro Peres-Neto

This repository contains code to reproduce all figures, tables and analyses included in the manuscript and in its supplementary information. This repository is directly linked to the Open Science Framework project for this study.

## Data availability

A stable data set containing all data required for analysis, reanalysis, and reproduction of figures and materials is available within the Dryad repository (accessible at https://doi.org/10.5061/dryad.rjdfn2zgj) related to this study.

## Keywords

biogeography, Chiroptera, climatic legacy, community assembly, paleoclimate, macroecology, phylogenetic relatedness, scale-dependence, diversification, spatial scales

## Abstract

<details>

<summary>Click to read</summary>

### Aim

Patterns of evolutionary relatedness among co-occurring species are driven by scale-dependent contemporary and historical processes. Yet, we still lack a detailed understanding of how these drivers impact the phylogenetic structure of biological communities. Here, we focused on bats -- one of the most speciose and vagile groups of mammals -- and test the predictions of three general biogeographical hypotheses that are particularly relevant to understanding how paleoclimatic stability, local diversification rates, and geographical scales shaped their present-day phylogenetic community structure.

### Location

Worldwide, across restrictive spatial extents: global, east-west hemispheres, biogeographical realms, tectonic plates, biomes, and ecoregions.

### Time period

Last Glacial Maximum (\~22,000 years ago) to present.

### Major taxa studied

Bats (Chiroptera).

### Methods

We estimated bat phylogenetic community structure across restrictive geographical extents and modelled it as a function of paleoclimatic stability, and *in situ* net diversification rates.

### Results

Limiting species pools from broader to local spatial scales Limiting geographical extents from larger to smaller scales strongly changed the phylogenetic structure of bat communities. The magnitude of these effects is less noticeable in the western hemisphere, where frequent among-realm biota interchange could have been maintained through bats adaptive traits. Highly phylogenetically related bat communities are generally more common in regions that changed less in climate since the last glacial maximum, supporting the expectation that stable climates allow for increased phylogenetic clustering. Finally, increased *in situ* net diversification rates are associated with greater phylogenetic clustering in bat communities.

### Main conclusions

We show that the worldwide phylogenetic structure of bat assemblages varies as a function of geographical extents, dispersal barriers, paleoclimatic stability and in situ diversification. The integrative framework used in our study, which can be applied to other taxonomic groups, has proven useful to not only explain the evolutionary dynamics of community assembly, but could also help tackle questions related to scale-dependence in community ecology and biogeography.

</details>

## Description of this research compendium

<details>

<summary><code> scripts/ </code></summary>

-   `S0_fun_BAMM.diversification.rate.estimation.R`: Function to configure BAMM diversification rate estimation files;
-   `S0_fun_bootstrap_logit_GLM_uncondition_quant.R`: Functions to compute and extract logistic regression coefficients;
-   `S0_fun_CommWeightedMeans.R`: Function to compute community weighted means;
-   `S0_fun_ggplot_theme_map.R`: Function to define `ggplot` theme for Figure 1;
-   `S0_fun_make_grid_sf.R`: Function to create a cell-grid over a shapefile of the world containing polygon layers;
-   `S0_fun_match_phylo_comm.R`: Function to match the species tips of a phylogenetic tree to a community matrix containing species occurrences across sites;
-   `S0_fun_ses.opt.rarefaction.phylostr.R`: Functions to compute the phylogenetic relatedness of communities using the traditional and a rarefaction approach based on the standardized effect sizes of mean phylogenetic pairwise and mean nearest taxon distances;
-   `S0_fun_ses.phylostr_non_parallelized_alternatives.R`: Function that implements expanded grids to calculate the standardized effect size for mean phylogenetic pairwise distances and mean nearest taxon distances (based on `PhyloMeasures::mpd.query()` and `PhyloMeasures::mntd.query()`);
-   `S0_fun_ses.phylostr.query.sf.R`: Functions that implements `PhyloMeasures::mpd.query()` and `PhyloMeasures::mntd.query()` within the framework of `picante::ses.mpd()` and `picante::ses.mntd()` to allow for faster parallel computation using SNOW/snowfall;
-   `S0_fun_sf.ses.phylostr.R`: Functions that modify `picante::ses.mpd()` and `picante::ses.mntd()` to allow for parallel computation using SNOW/snowfall;
-   `S0_HypothesesRepresentationFigures.R`: Routine to simulate and generate the figures that represent the hypotheses being tested and that are inserted within Table 1;
-   `S00_REnvironmentPreparation.R`: Routine to prepare (install and load packages) the R environment for analysis and reanalysis of the data in this study;
-   `S1a_50KM_DatasetPreparation.R`: Routine to prepare the geographical data set for this study;
-   `S1b_Chiroptera_Comm_Phylo_DatasetPreparation.R`: Routine to prepare the community presence absence data, the phylogenetic relationship hypothesis and the maximum crade credibility tree used in our study;
-   `S1c_Climate_Contemp_LGM_DatasetPreparation.R`: Routine to prepare the contemporary climatic, paleoclimatic data, and climatic legacies;
-   `S2a_PhyloStr_SamplingPool_MPD_MNTD_mod.ses.mpd.query.sf.R`: Routine to apply a null-model framework to compute the phylogenetic structure of communities across a gradient of restrictive geographical extents;\
-   `S2b_Summary_Statistics_NRI_NTI_Table_S2.1.R`: Code to summarise statistics for the phylogenetic structure of communities, and to create Table S2.1;
-   `S2c_mergePhyloStructure_rmanova_rmmcp_Table_S2.2.R`: Code to merge data on community phylogenetic structure, perform robust repeated measurement analyses of variance, and create Table S2.2.;
-   `S2d_PhyloStr_SamplingPool_MPD_MNTD_P_Z_Combining_Figure_2.R`: Code to perform Stouffer's meta-analytic probability combination tests on the probabilities from the computations of indices for community phylogenetic relatedness;
-   `S2e_PhyloStr_SamplingPool_MPD_MNTD_PhyloSamples.R`: Routine to apply a null-model framework to compute the phylogenetic structure of communities across a gradient of restrictive geographical extents across phylogenetic trees randomly sampled from the phylogenetic relationship hypothesis used in this study;
-   `S2f_PhyloStr_SamplingPool_MPD_MNTD_Fixed_Rarefaction.R`: Routine to perform a rarefaction-based adjustment for biases introduced by differences in sizes of species richness under the null-model framework to compute the phylogenetic structure of communities across a gradient of restrictive geographical extents across phylogenetic trees. This code allows for the selection of fixed species richness across all communities;
-   `S2g_PhyloStr_SamplingPool_MPD_MNTD_Relative_Rarefaction.R`: Routine to perform a rarefaction-based adjustment for biases introduced by differences in sizes of species richness under the null-model framework to compute the phylogenetic structure of communities across a gradient of restrictive geographical extents across phylogenetic trees. This routine allows for the selection of relative species richness across all communities;
-   `S2h_PhyloStr_SamplingPool_MPD_MNTD_Relative_Min_Size_Rarefaction.R`: Routine to perform a rarefaction-based adjustment for biases introduced by differences in sizes of species richness under the null-model framework to compute the phylogenetic structure of communities across a gradient of restrictive geographical extents across phylogenetic trees. Here, the bias is adjusted by repeatedly randomly subsampling (rarefying) any given local community matrix to have the same number of species as the immediately inferior nested geographical extent;
-   `S2i_PhyloStr_SamplingPool_Representation_Map_Figure_1.R`: Code to represent community phylogenetic structure calculated for several geographical scales into maps, creating Figure 1;
-   `S3a_BAMM_netDiv_rate_estimation.R`: Code to generate the control file, assess and estimate net diversification rates for each species using Bayesian Analyses for Macroevolutionary Mixtures using the maximum crade credibility tree computed from the phylogenetic hypothesis used here;
-   `S3b_BAMM_netDiv_rate_estimation_PhyloSamples.R`: Code to generate the control file, assess and estimate net diversification rates for each species using Bayesian Analyses for Macroevolutionary Mixtures for phylogenetic trees randomly sampled from the phylogenetic hypothesis used in this study;
-   `S3c_samplingProbabilities_Table_S2.3.R`: Code to compute sampling probabilities within the phylogenetic hypothesis used in this study, and to create Table S2.3;
-   `S4a_CWM_netDiv.R`: Code to compute weighted means for diversification rates across all communities;
-   `S4b_CWM_netDiv_PhyloSamples.R`: Code to compute average weighted means for diversification rates across all communities that were obtained for each phylogenetic tree randomly sampled for the phylogenetic hypothesis used in this study;
-   `S5_PhyloStr_Descriptive_Quantiles_Figures_3_4.R`: Code to plot the mean phylogenetic relatedness of bat communities (for both NRI and NTI) across each percentile of the predictor variables of interest (i.e., historical change in temperature, historical change in precipitation and in situ net diversification rates), and generate Figures 3 and 4;
-   `S6_PhyloStr_Logistic_Bootstrap_Quantiles_Figure_5_Table_2.R`: Code to applu conditionally unbiased bounded influence robust logistic regressions to test how changes in historical climatic stability and in situ diversification rates independently increased (or decreased) the likelihood of a community being composed of highly phylogenetically related species. This code creates Figure 5 and Table 2;
-   `S7a_CHELSA_TraCE21k_download.R`: Code to download and decompress climatic data from CHELSA covering the current period until 22,000 years ago, in intervals of 500 years;
-   `S7b_CHELSA_TraCE21k_VoCC_extract.R`: Code to calculate gradient-based change velocities in temperature and in precipitation from the LGM to the contemporary period and assessed how they influenced the phylogenetic relatedness of bat communities across geographical scales.

</details>

<details>

<summary><code> manuscript/supp_info/ </code></summary>

</details>

<details>

<summary><code> data/ </code></summary>

</details>

*Will be included soon.*

## Citing

Study:
