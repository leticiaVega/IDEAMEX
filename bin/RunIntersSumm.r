#
#
#
#
library(limma)

IntersectionSummary<-function(fnIntersectionTable)
{
    fnColSummay<-as.integer(colSums(fnIntersectionTable))
    fnColSummay<-c(fnColSummay,sum(fnColSummay))
    fnAllTableSumm<-vennCount(fnIntersectionTable)
    fnAllTableSumm<-fnAllTableSumm[seq(nrow(fnAllTableSumm),1),]
    fnAllTableSumm["1",]<-fnColSummay
    return(fnAllTableSumm)
}