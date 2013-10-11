##Programmers: Daniel Bottomly with assistance from Peter Ryabinin
##A set of functions to prioritize genes based on one or more samples being differentially expressed

default.weight.func=function(x) as.matrix(dist(t(x)))

setGeneric("outlyingDegree", def=function(obj,...) standardGeneric("outlyingDegree"))
setMethod("outlyingDegree", signature("ExpressionSet"), function(obj, k, type=c("non.weight", "weight.before", "weight.after"), weight.func=default.weight.func)
          {
                type = match.arg(type)
                
                if(!is.function(weight.func))
                {
                    stop("ERROR: weight.func needs to be a function")
                }
                
                use.exprs <- exprs(obj)
                add.params <- list(od.k=k, weight.func=weight.func)
                
                return(switch(type, non.weight=od.matrix(use.exprs, list(od.k=k)), weight.before=weight.od.b(use.exprs, add.params), weight.after=weight.od.a(use.exprs, add.params)))
          })

#these two based off of the original POD and spatial weighted outlier detection by Kou, Lu and Chen
weight.od.b <- function(exprs, add.params)
{
    check.params(add.params, c(od.k="numeric", weight.func="function"))
    check.mat(exprs)
    
    k <- add.params$od.k
    weight.func <- add.params$weight.func
    
    samp.weights <- weight.func(exprs)
    
    check.mat(samp.weights)
    
    samp.weights <- samp.weights[colnames(exprs), colnames(exprs)]
    
    if (k < 1 || k >= ncol(exprs))
    {
        stop("ERROR: od.k should be between 1 and ncol(exprs)-1")
    }
    
    pat.mat <- sapply(1:ncol(exprs), function(x)
       {
            sub.samp.weights <- samp.weights[x,-x]/sum(samp.weights[x,-x])
            temp.dist <- abs(exprs[,x] - exprs[,-x])*matrix(sub.samp.weights, nrow=nrow(exprs), ncol=length(sub.samp.weights), byrow=TRUE)

            return(apply(temp.dist, 1, function(y) sum(sort(y, decreasing=FALSE)[1:k])))
        
      })
          
    colnames(pat.mat) <- colnames(exprs)
    
    return(pat.mat)
}

weight.od.a <- function(exprs, add.params)
{
    check.params(add.params, c(od.k="numeric", weight.func="function"))
    check.mat(exprs)
    
    k <- add.params$od.k
    weight.func <- add.params$weight.func
    
    samp.weights <- weight.func(exprs)
    
    check.mat(samp.weights)
    
    samp.weights <- samp.weights[colnames(exprs), colnames(exprs)]
    
    if (k < 1 || k >= ncol(exprs))
    {
        stop("ERROR: od.k should be between 1 and ncol(exprs)-1")
    }
    
    pat.mat <- sapply(1:ncol(exprs), function(x)
       {
            temp.dist <- abs(exprs[,x] - exprs[,-x])

            return(apply(temp.dist, 1, function(y)
                         {
                            ord.y <- order(y, decreasing=FALSE)[1:k]
                            return(sum(y[ord.y]*(samp.weights[x,-x][ord.y]/sum(samp.weights[x,-x][ord.y]))))
                         }))
        
      })
          
    colnames(pat.mat) <- colnames(exprs)
    
    return(pat.mat)
}

#where check.names is a named character vector of the form: name=type (e.g. od.k="integer")
check.params <- function(param.list, check.names)
{
    if ((class(param.list) == "list" && all(names(check.names) %in% names(param.list))) == FALSE)
    {
        stop(paste0("ERROR: param.list needs to be a named list with elements:", paste(names(check.names), collapse=",")))
    }
    
    for (i in names(check.names))
    {
        if (class(param.list[[i]]) != check.names[i])
        {
            stop(paste("ERROR: param.list element-type mismatch:", i, "is not of type", check.names[i]))
        }
    }
}

check.mat <- function(mat)
{
    if (is.null(dimnames(mat)) || is.null(dimnames(mat)[[1]]) || is.null(dimnames(mat)[[2]]))
    {
        stop("ERROR: all matrices must have row and column names")
    }
}

od.matrix <- function(exprs, add.params)
{
    
    hollow.mat.func <- function(x)
    {
        hollow.mat <- 1-diag(nrow=ncol(x))
        dimnames(hollow.mat) <- list(colnames(x), colnames(x))
        return(hollow.mat)
    }
    
    add.params$weight.func <- hollow.mat.func
    
    scaled.od.vals <- weight.od.b(exprs, add.params)
    
    return(scaled.od.vals*(ncol(exprs)-1))
    
}

#Simple Z score computation, keep add.params set to NULL
zscore.matrix <- function(inp.mat, add.params)
{
        mat.z <- t(scale(t(inp.mat), center=TRUE, scale=TRUE))
        return(mat.z)
}

robust.z <- function(exprs, add.params)
{
    #as defined in tibshurani and hastie 2007 biostats
    
    base.mad <- apply(exprs, 1, mad)

    ygi <- t(apply(exprs, 1, function(x) (x-median(x))))

    ygi.cor <- ygi/base.mad
    
    return(ygi.cor)
}

#Compute the P-values as in Emerson and Emerson 2011
get.sim.p.half.tie <- function(rank.vec, control.pos)
{
    abs.rank.vec <- abs(rank.vec)
    return((sum(abs.rank.vec[-control.pos] > abs.rank.vec[control.pos]) + ((1/2)*(sum(abs.rank.vec[-control.pos] == abs.rank.vec[control.pos]))))/(length(rank.vec)-1))
}