\name{Methods for patient-specific outlier detection}
\alias{outlyingDegree}
\alias{outlyingDegree,matrix-method}
\alias{outlyingDegree,ExpressionSet-method}
\alias{zScore}
\alias{zScore,matrix-method}
\alias{zScore,ExpressionSet-method}
\title{
	Methods for patient-specific outlier detection
}
\description{
	Implementations of the outlying degree and Zscore methodologies desribed in Bottomly et al 2013.
}
\usage{
\S4method{outlyingDegree}{ExpressionSet}(obj, k, type=c("non.weight", "weight.before", "weight.after"), weight.func=default.weight.func)
\S4method{outlyingDegree}{matrix}(obj, k, type=c("non.weight", "weight.before", "weight.after"), weight.func=default.weight.func) 
\S4method{zScore}{ExpressionSet}(obj, robust=FALSE) 
\S4method{zScore}{matrix}(obj, robust=FALSE)
}
\arguments{
	\item{obj}{Either a \code{matrix} with rows corresponding to genes and columns corresponding to samples or an \code{ExpressionSet} object }
	\item{k}{An integer value that can range from 1 to ncol(obj)-1}
	\item{type}{One of either "non.weight", "weight.before" or "weight.after" which corresponds to the standard outlying degree method, whether to apply the weights specified in \code{weight.func} before computing the distance or after computing the distance respectively.}
	\item{weight.func}{A function returning a square matrix containing the weights corresponding to the sample-sample dissimilarities.}
	\item{robust}{Whether a robust Zscore should be computed e.g. (x-median(x))/mad(x)}
	}
	
\author{
Daniel Bottomly
}

\examples{
data(simEset)
out.deg.eset <- outlyingDegree(cur.eset, k=1)
cur.mat <- exprs(simEset)
out.deg.mat <- outlyingDegree(cur.mat, k=1)
all.equal(out.deg.eset, out.deg.mat)

zscore.eset <- zScore(cur.eset)
zscore.mat <- zScore(cur.mat)

all.equal(zscore.eset, zscore.mat)
}
	
\keyword{utilities}
