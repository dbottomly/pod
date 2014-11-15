set.seed(123)
require(Biobase)

simExprs <- ExpressionSet(assayData=matrix(rnorm(1000, mean=6, sd=2), nrow=100, ncol=10, dimnames=list(paste("gene", 1:1000,sep="_"),paste("sample", 1:10 ,sep="_"))))
