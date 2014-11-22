set.seed(123)
require(Biobase)

simEset <- ExpressionSet(assayData=matrix(rnorm(10000, mean=6, sd=2), nrow=1000, ncol=10, dimnames=list(paste("gene", 1:1000,sep="_"),paste("sample", 1:10 ,sep="_"))))
