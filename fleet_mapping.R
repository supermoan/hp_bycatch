# @Author André Moan
# @Email andre.moan@hi.no
# @Description This script contains the part of a larger analysis that is concerned 
#              with mapping ordinary (non-reference) fishing vessels to reference 
#              vessels. The purpose is to make bycatch predictions using generalized
#              additive models (GAMs) that include random effect terms possible. 
#
#              The purpose of this R script is to document the process by which
#              we included random effects in predictions for GAMs.
#
#              Note that vessel names have been redacted for privacy purposes.

# We start out with a data frame, x, that contains fishing effort data for all vessels,
# (both ordinary and reference vessels) aggregated by region and fishery. The contents of
# x that is shown below reflect the data that was used in producing bycatch estimates 
# given in Moan, Skern-Mauritzen, Vølstad & Bjørge. New estimates of bycatch of harbour 
# porpoise (Phocoena phocoena) in Norwegian gillnet fisheries suggest unsustainable 
# incidental mortality (in submission). Other data aggregation schemes were also
# attempted. 

# > str(x)
#'data.frame':	19868 obs. of  9 variables:
# $ region    : num  3 3 3 3 3 2 3 4 2 3 ...
# $ fishery   : Factor w/ 3 levels "other","angler",..: 1 2 3 1 2 3 3 2 1 1 ...
# $ reference : num  0 0 0 0 0 0 0 0 0 0 ...
# $ hauls     : num  5 1 5 27 186 17 191 1 1 16 ...
# $ catch     : num  570.6 39.6 1886.2 9316.9 28027.7 ...
# $ vesselsize: num  6.8 6.8 6.8 10.6 10.6 ...
# $ vessel    : chr  "REDACTED1" "REDACTED1" "REDACTED1" "REDACTED2" ...
# $ effort    : num  5 1 5 27 186 17 191 1 1 16 ...

# get_ranks(): takes a data frame containing the columns vessel, vesselsize and effort. Additionally,
# each variable specified in vars must also be a column in data. The data frame is first split by
# vessel, and then each subset is ordered decrementally by effort. Labels for each effort datum
# is generated from combinations of the stratifying variables specified. Returns a list of lists,
# where the inner list is comprised of three vectors, vessel, size and effort. 
get_ranks <- function(data, vars = c("fishery", "region")) {
    lapply(split(data, data$vessel), function(x) {
        x <- x[order(x$effort, decreasing=T),]
        labels <- do.call(paste0, x[,vars])
        list(vessel = x$vessel[1], 
             size = x$vesselsize[1],
             effort = setNames(x$effort, labels)
        )
    })
}

ranks.ref <- get_ranks(x[x$reference == 1,]) # ranks for the reference vessels
ranks.all <- get_ranks(x[x$reference == 0,]) # ranks for the non-reference vessels

# > str(head(ranks.ref, n = 3))
#List of 3
# $ REDACTED.ref1:List of 3
# ..$ vessel: chr "REDACTED.ref1"
# ..$ size  : num 12.2
# ..$ effort: Named num [1:3] 204 44 10
# .. ..- attr(*, "names")= chr [1:3] "angler4" "other4" "cod4"
# $ REDACTED.ref2:List of 3
# ..$ vessel: chr "REDACTED.ref1"
# ..$ size  : num 14.7
# ..$ effort: Named num [1:3] 215 132 32
# .. ..- attr(*, "names")= chr [1:3] "cod1" "angler1" "other1"
# $ REDACTED.ref3:List of 3
# ..$ vessel: chr "REDACTED.ref3"
# ..$ size  : num 12.8
# ..$ effort: Named num [1:4] 72 69 17 2
# .. ..- attr(*, "names")= chr [1:4] "other3" "cod2" "other4" "cod4"

# > str(head(ranks.all, n = 3))
#List of 3
# $ REDACTED1:List of 3
# ..$ vessel: chr "REDACTED1"
# ..$ size  : num 6.8
# ..$ effort: Named num [1:3] 5 5 1
# .. ..- attr(*, "names")= chr [1:3] "other3" "cod3" "angler3"
# $ REDACTED2:List of 3
# ..$ vessel: chr "REDACTED2"
# ..$ size  : num 10.6
# ..$ effort: Named num [1:4] 191 186 27 17
# .. ..- attr(*, "names")= chr [1:4] "cod3" "angler3" "other3" "cod2"
# $ REDACTED3:List of 3
# ..$ vessel: chr "REDACTED3"
# ..$ size  : num 10.6
# ..$ effort: Named num [1:3] 227 61 7
# .. ..- attr(*, "names")= chr [1:3] "cod1" "other1" "angler1"

# iterate over the list of non-reference vessels and compare data for each non-reference vessel
# with data for each reference vessels. In each iteration, return exactly one reference vessel 
# name. 
m <- sapply(ranks.all, function(x) {
    
    # determine eligible reference vessels
    s <- sapply(ranks.ref, function(x) x$size)
    for (threshold in seq(2,10,1)) {
        i <- which(abs(s-x$size) < threshold)
        if (length(i)) break
    }
    
    # find the best match amongst eligible vessels
    scores <- sapply(ranks.ref[i], function(y) rbo(x$effort, y$effort, p = 0.8))
    
    # in case of a tie, pick a random vessel amongst tied vessels
    sample(names(scores)[all.max(scores)], size = 1)
})

# > str(m)
# Named chr [1:8188] "DEDACTED.ref5" "DEDACTED.ref5" "DEDACTED.ref10" "Vågøybuen" "DEDACTED.ref4" "DEDACTED.ref4" "DEDACTED.ref3" "DEDACTED.ref11" ...
# - attr(*, "names")= chr [1:8188] "DEDACTED1" "DEDACTED2" "DEDACTED3" "DEDACTED4" ...

# overwrite original non-reference vessel names in newdata with most closely matching reference vessels
newdata$fvessel <- factor(m[match(newdata$vessel, names(m))])

# now we can pass newdata to predict.gam to predict with random effects!
