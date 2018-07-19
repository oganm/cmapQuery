library(memoise)
library(ConnectivityMap)
load_all()

d = 100000

data("instances")
data("rankMatrix")
# when memoising data types are important
set.seed(1)
preCalcKs = preCalcRandomKs(instances$cmap_name,ncol(rankMatrix))

saveRDS(preCalcKs,'data-raw/legacyPreCalc.rds')
