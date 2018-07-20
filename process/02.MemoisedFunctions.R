library(memoise)
library(ConnectivityMap)
load_all()

d = 100000

data("instances")
data("rankMatrix")
# when memoising data types are important
set.seed(1)
cmapPreCalc = preCalcRandomKs(instances$cmap_name,ncol(rankMatrix))

devtools::use_data(cmapPreCalc)

# saveRDS(preCalcKs,'data-raw/legacyPreCalc.rds')



legacyMsigDBPreCalc = specificityPreCalculation(signatures = MSigDBLegacy,
                                                rankMatrix = rankMatrix,
                                                chems = instances$cmap_name,
                                                preCalc = cmapPreCalc,
                                                cores = 5)

devtools::use_data(legacyMsigDBPreCalc)
