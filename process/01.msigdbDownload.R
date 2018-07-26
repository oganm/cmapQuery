library(magrittr)
library(purrr)
library(stringr)

# download original msigDB data from broad institute ------
download.file('https://portals.broadinstitute.org/cmap/msigdb_gene_sets.zip',destfile = 'data-raw/msigdb_gene_sets.zip')
unzip('data-raw/msigdb_gene_sets.zip',exdir = 'data-raw/')
MsigUp = readLines('data-raw/msigdb_up_mapped_to_HG_U133A.gmt') %>% strsplit('\t')
names = MsigUp %>% map_chr(1) %>% gsub(pattern = '_UP',replacement = '',x=.)

MsigDown = readLines('data-raw/msigdb_dn_mapped_to_HG_U133A.gmt')%>% strsplit('\t')

names2 = MsigDown %>% map_chr(1) %>% gsub(pattern = '_DN',replacement = '',x=.)

assertthat::see_if(all(names==names2))

MSigDBLegacy = lapply(1:length(MsigUp), function(i){
    upTags = MsigUp[[i]][c(-1,-2)]
    downTags = MsigDown[[i]][c(-1,-2)]
    return(list(upTags = upTags,downTags = downTags))
})
names(MSigDBLegacy) = names
devtools::use_data(MSigDBLegacy,overwrite = TRUE)

# new msigdb data -------------
library(MSigDB) # github.com/oganm/MsigDB
# find up-down tagged lists

sigNames = MSigDB %>% unlist(recursive = FALSE) %>% names

upSigs = sigNames[sigNames %>% grep('UP$',.,ignore.case = TRUE)] %>% str_extract('.*?(?=_UP$)')
downSigs = sigNames[sigNames %>% grep('DN$',.,ignore.case = TRUE)] %>% str_extract('.*?(?=_DN$)')

upDown = intersect(upSigs,downSigs)

MSigDB62 = lapply(str_split(upDown,'\\.'),function(x){
    upTags = MSigDB[[x[1]]][[paste0(x[2],'_UP')]]
    downTags = MSigDB[[x[1]]][[paste0(x[2],'_DN')]]
    return(list(upTags = upTags, downTags = downTags))
})
names(MSigDB62) = upDown

devtools::use_data(MSigDB62,overwrite = TRUE)
