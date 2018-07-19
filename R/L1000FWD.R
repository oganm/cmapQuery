# currently broken
L1000FWD = function(upGenes,downGenes){
    L1000FWD_URL <- 'http://amp.pharm.mssm.edu/L1000FWD/'

    payload = list(up_genes = upGenes,
                   down_genes = downGenes)
    response <- httr::POST(paste0(L1000FWD_URL, 'sig_search/'), body=payload, encode='json')

}
