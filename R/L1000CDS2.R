#' @export
L1000CDS2 = function(upGenes= NULL, downGenes = NULL, signature = NULL,
                     mimic = TRUE, combination = TRUE,dbVersion = 'latest'){
    baseLink = 'http://amp.pharm.mssm.edu/L1000CDS2/query'
    if(is.null(signature)){
        assertthat::assert_that(!is.null(upGenes) & !is.null(downGenes))
        payload = list(data = list(upGenes = upGenes,
                                   dnGenes = downGenes),
                       config = list(aggravate = mimic,
                                     searchMethod = 'geneSet',
                                     share = FALSE,
                                     combination = combination,
                                     `db-version` = dbVersion))

        response = httr::POST(baseLink,body = payload,encode = 'json')
    } else{
        payload = list(data  = list(genes = names(signature),vals = unlist(unlist(signature))),
                       config = list(aggravate = mimic,
                                     searchMethod = 'CD',
                                     share = FALSE,
                                     combination = combination,
                                     `db-version` = dbVersion))
        response = httr::POST(baseLink,body = payload,encode = 'json')

    }

    response <- jsonlite::fromJSON(httr::content(response, 'text'))
    return(response)
}
