#' @export
a =function(V,n){
    V %<>% as.matrix
    t = nrow(V)
    if(nrow(V)==0){
        return(0)
    }
    ((1:t)/t - V/n) %>% matrixStats::colMaxs()
}


#' @export
b = function(V,n){
    V %<>% as.matrix
    t = nrow(V)
    if(nrow(V)==0){
        return(0)
    }
    (V/n - (1:t - 1)/t) %>% as.matrix %>% matrixStats::colMaxs()
}


#' @export
ks = function(a,b){
    (a>b)*a + (b>a)*(-b)
}

#' @export
ksCalc = function(V,n){
    a = a(V,n)
    b = b(V,n)
    return(ks(a,b))
}

#' #' @export
#' memoKsCalc = memoise::memoise(ksCalc)


#' @export
randomV = function(length,n,d){
    replicate(d,sample(x=1:n,
                       size = length,replace = FALSE) %>% sort)
}

#' #' @export
#' memoRandomV = memoise::memoise(randomV)
#'

#' @export
directRandomKsCalc = function(length,n,d){
    randomV = randomV(length,n,d)
    ksCalc(randomV,n)
}

#'
#' #' @export
#' memoDirectRandomKsCalc = memoise::memoise(directRandomKsCalc)


#' @export
score = function(kUp,kDown){
    (sign(kUp)!=sign(kDown))*(kUp-kDown)
}


#' @export
scoreCalc = function(Vup,Vdown,n){
    aUp = a(Vup,n)
    bUp = b(Vup,n)
    kUp = ks(aUp,bUp)

    aDown = a(Vdown,n)
    bDown = b(Vdown,n)
    kDown = ks(aDown,bDown)
    score = score(kUp,kDown)
    return(data.frame(kUp = kUp, kDown = kDown,score = score,instance= colnames(Vup)))
}

#' @export
connectivityMapEnrichment = function(upTags,downTags,rankMatrix,chems,pAdjustMethod = 'fdr' ,d=100000,
                                     preCalc = NULL){

    if(!is.null(preCalc)){
        memoDirectRandomKsCalc = preCalc
    }
    assertthat::assert_that(ncol(rankMatrix)==length(chems))

    n = rankMatrix %>% nrow
    Vup = rankMatrix[upTags,] %>% apply(2,sort)
    Vdown = rankMatrix[downTags,] %>% apply(2,sort)
    scores = scoreCalc(Vup,Vdown,n)

    # "p to be max( si) and q to be min( si) across all instances in the collection c"
    p = max(scores$score)
    q = min(scores$score)

    scores %<>% dplyr::mutate(ConScore = (score>0)*(score/p) + (score<0)*(-score/q))

    # "The Kolmogorov-Smirnov statistic is computed for the set of t instances
    # in the list of all n instances in a result ordered in descending order of
    # connectivity score and up score (see how connectivity score is
    # calculated), giving an enrichment score ks0."
    scores = scores %>%
        dplyr::mutate(order = seq_len(nrow(scores))) %>%
        dplyr::arrange(desc(ConScore),desc(kUp))

    uniqueChems = chems %>% unique

    confidence = uniqueChems %>% sapply(function(chem){
        # browser()
        # print(chem)
        chemInstances = which(chems %in% chem)

        relevantInstances = scores %>% dplyr::filter(order %in% chemInstances)

        relevantUpCount = (relevantInstances$score>0) %>% sum
        relevantDownCount = (relevantInstances$score<0) %>% sum
        relevantMehCount = (relevantInstances$score==0) %>% sum
        nonNull = (relevantUpCount>=relevantDownCount)*relevantUpCount +
            (relevantUpCount<relevantDownCount)*(relevantDownCount)
        nonNull = nonNull/nrow(relevantInstances)

        # data.table::fsort
        V = match(chemInstances,scores$order) %>% sort
        ks0 = ksCalc(V,ncol(rankMatrix))

        # if precalculated permutations are entered, use those
        if(!is.null(preCalc)){
            ksPerm = memoDirectRandomKsCalc(length(chemInstances),ncol(rankMatrix),d)
            # Vrandoms = memoRandomV(length(chemInstances),ncol(rankMatrix),d)
            # ksPerm = memoKsCalc(Vrandoms,ncol(rankMatrix))
        } else{
            Vrandoms = randomV(length(chemInstances),ncol(rankMatrix),d)
            ksPerm = ksCalc(Vrandoms,ncol(rankMatrix))
        }


        q = sum(abs(ksPerm) >= abs(ks0))
        p = q/d

        return(c(enrichment = ks0,
                 p = p,
                 nonNull = nonNull,
                 instanceCount =  length(chemInstances)))
    }) %>% t

    confidence %<>% as.data.frame
    confidence$FDR =  p.adjust(confidence$p,method = pAdjustMethod)
    confidence = confidence[c('enrichment','p','FDR','nonNull','instanceCount')]

    return(list(instanceScores = as.data.frame(scores), chemScores = confidence))
}


# precalculates random Ks. chems is the vector of unique chemicals, experiment count is the number
# of experiments (ncol(rankMatrix)) d is the permuation count.
# the output is a memoised function that should be given to connectivityMapEnrichment
# to speed things up
#' @export
preCalcRandomKs = function(chems,experimentCount, d = 100000,preCalc = NULL, debug = FALSE){
    if(!is.null(preCalc)){
        memoDirectRandomKsCalc = preCalc
    } else{
        memoDirectRandomKsCalc = memoise::memoise(directRandomKsCalc)
    }
    instanceLengths = table(chems) %>% unique

    if(!debug){
        pb = txtProgressBar(min = 0, max = length(instanceLengths), initial = 0)
    }
    iLens = seq_along(instanceLengths)
    alkKs = sapply(iLens, function(i){
        l = instanceLengths[i]
        memoDirectRandomKsCalc(l,experimentCount,d)

        if(debug){
            print(glue::glue('({l},{experimentCount},{d})'))
            print(glue::glue('({class(l)},{class(experimentCount)},{class(d)})'))
        } else{
            setTxtProgressBar(pb,i)
        }
    })
    return(memoDirectRandomKsCalc)
}
