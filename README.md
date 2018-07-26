
Readme
------

This repository aims to replicate the analysis performed at Lamb et al 2006. but allow using the L1000 data instead of the original Cmap data.

This is currently more of a "notes to self" rather than a readme.

Working example
---------------

``` r
devtools::load_all()
```

    ## Loading cmapQuery

``` r
library(ConnectivityMap) # bioconductor houses the old cmap data
data("rankMatrix")
data("instances")

library(parallel)
library(gemmaAPI) # devtools::install_github('PavlidisLab/GemmaAPI.R') used for gene-platform matching

# the data in the connectivity map is marked with probes so our input has to be 
# probes to. The platform that is used is GPL96. New data comes with gene symbols
# so these steps won't be necesarry 

upGenes = c("BFSP1", "C11orf63", "LCN2", "SIX2", "SIX2", "TLR2")

downGenes = c("MUC5B", "MUC5B", "TECTA", "COLQ", "FAM124B", "ASB4")

gpl96 = getAnnotation('GPL96')

upGenes %<>% annotationProbesetMatch(gpl96,removeNAs = TRUE)
downGenes %<>% annotationProbesetMatch(gpl96, removeNAs = TRUE)

# now pre-train to make everything go faster later on
# this package already includes precalculated data for
# the old cmap data. not running here because it takes a while
# considering using a smaller d than the default one to get things going
# set.seed(1)
# cmapPreCalc = preCalcRandomKs(instances$cmap_name,ncol(rankMatrix))

analysis = connectivityMapEnrichment(upGenes,
                                     downGenes,
                                     rankMatrix,
                                     instances$cmap_name,
                                     preCalc = cmapPreCalc)

# see the scores given for individual chemicals
analysis$chemScores %>% head()
```

    ##                   enrichment       p          FDR   nonNull instanceCount
    ## metformin        -0.16524590 0.90742 0.9849193864 0.2000000            10
    ## phenformin        0.37386417 0.21757 0.7933123398 0.1428571             7
    ## phenyl biguanide -0.60672131 0.00001 0.0002727083 0.0000000             1
    ## valproic acid     0.15912568 0.09739 0.6629342784 0.2631579            57
    ## estradiol         0.07991582 0.95487 0.9928033229 0.2432432            37
    ## alpha-estradiol   0.20696721 0.43703 0.8859794752 0.2500000            16
    ##                  reliable
    ## metformin           FALSE
    ## phenformin          FALSE
    ## phenyl biguanide    FALSE
    ## valproic acid       FALSE
    ## estradiol           FALSE
    ## alpha-estradiol     FALSE

``` r
#the analysis also includes reliability and specificity scores. this
#historically uses a frozen version of MsigDB database. Both legacy MsigDB
#database used for CMAP and a new list based on current MSigDB is included
#within the package (MSigDBLegacy, MsigDB62). Legacy cmap comes as probe names,
#new list comes as gene names

# getting specificity and reliability requires pre-training with the MSigDB data
# this takes some time so it can be parallelized by the cores argument
# legacyMsigDBPreCalc = specificityPreCalculation(signatures = MSigDBLegacy,
#                                                 rankMatrix = rankMatrix,
#                                                 chems = instances$cmap_name,
#                                                 preCalc = cmapPreCalc,
#                                                 cores = 5)

# this package includes legacyMsigDBPreCalc now. Adding the precalcs for the new data
# will not be feasible though. recommended option is to load them as needed from an rds
# experiment with compression options



analysis$chemScores$specificity = specificityCalc(analysis$chemScores, legacyMsigDBPreCalc)

analysis$chemScores %>% head
```

    ##                   enrichment       p          FDR   nonNull instanceCount
    ## metformin        -0.16524590 0.90742 0.9849193864 0.2000000            10
    ## phenformin        0.37386417 0.21757 0.7933123398 0.1428571             7
    ## phenyl biguanide -0.60672131 0.00001 0.0002727083 0.0000000             1
    ## valproic acid     0.15912568 0.09739 0.6629342784 0.2631579            57
    ## estradiol         0.07991582 0.95487 0.9928033229 0.2432432            37
    ## alpha-estradiol   0.20696721 0.43703 0.8859794752 0.2500000            16
    ##                  reliable specificity
    ## metformin           FALSE   0.9205298
    ## phenformin          FALSE   0.1860465
    ## phenyl biguanide    FALSE   0.8389831
    ## valproic acid       FALSE   0.5328947
    ## estradiol           FALSE   1.0000000
    ## alpha-estradiol     FALSE   0.7364341

Progress / Notes to self
------------------------

This repository aims to replicate the analysis performed at Lamb et al 2006. but allow using the L1000 data instead of the original Cmap data. With L1000, the data is much larger and there are a few different version of the stage 5 data around created by different approaches

-   The real stage 5 data by LINCS people
-   Our stage 5 data calculated by comparing to controls
-   Avi Mayan's stage 5 data, calculated by comparing to controls from the same lane. it's a little confusing. comunications pending

Since I want to experiment with all of those, this needs to be dataset independent. However my original methodology relied on memoisation for speed which is problematic if the background data can change. I need to decide when to do this memoisation and how much should I expose in on user side

Current process is to run `preCalcRandomKs` before running the experiment manually to get a memoised function and use it as an input to `connectivityMapEnrichment`. This will take a long time for the new L1000 data. I can start to see why Mayan's group abondoned this method.

To be able to quickly compare the studies it might make sense to use Mayan's web apps. I am implementing wrappers for the API calls of the apps.
