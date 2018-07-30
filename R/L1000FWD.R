#' @export
L1000FWD = function(upGenes,downGenes){
    L1000FWD_URL <- 'http://amp.pharm.mssm.edu/L1000FWD/'

    payload = list(up_genes = upGenes,
                   down_genes = downGenes)
    response <- httr::POST(paste0(L1000FWD_URL, 'sig_search'), body=payload, encode='json')

    response <- jsonlite::fromJSON(httr::content(response, 'text'))

    response = httr::GET(paste0(L1000FWD_URL,'result/topn/',response$result_id))

    response <- jsonlite::fromJSON(httr::content(response, 'text'))



    response %<>% lapply(function(x){
        chemInfo = x$sig_id %>% sapply(function(y){
            response = httr::GET(paste0(L1000FWD_URL,'sig/',y))
            response <- jsonlite::fromJSON(httr::content(response, 'text') %>% gsub('NaN','"NA"',.))

            return(c(pert_desc = response$pert_desc,
                     pert_dose = response$pert_dose,
                     cell_id = response$cell_id,
                     pert_id = response$pert_id,
                     pert_time = response$pert_time
                     ))
        })

        cbind(x, t(chemInfo))
    })

}

