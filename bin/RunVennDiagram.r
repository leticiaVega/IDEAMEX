# tabla<-read.table("~/Temp/Ejercicio2.txt",sep="\t",quote="",header=T,row.names=1)
# Shrubland<-subset(tabla,Shrubland>0,5)
# ......
# fnData<-list(row.names(Untreated_WW),row.names(Treated_WW),row.names(Freshwater),row.names(Rainfed),row.names(Shrubland))
# names(fnData)<-c("Untreated_WW","Treated_WW","Freshwater","Rainfed","Shrubland")
# source("~/bin/DifferentialExp/RunVennDiagram.r")
# VennDiag(fnData,"~/Temp","Ejercicio")

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
if(!("package:VennDiagram" %in% search()))
{
    library(VennDiagram)
}
if(!("package:limma" %in% search()))
{
    library(limma)
}
if(!("package:UpSetR" %in% search()))
{
    library(UpSetR)
}

if(!("package:ggplot2" %in% search()))
{
    library(ggplot2)
}


IntersectionSummary<-function(fnIntersectionTable,fnFileName)
{
    fnAllTableSumm<-vennCounts(fnIntersectionTable)
    fnAllTableSumm<-fnAllTableSumm[seq(nrow(fnAllTableSumm),1),]
    fnColSumm<-as.integer(colSums(fnIntersectionTable))
    fnWeightMatrix<-matrix(0,nrow=nrow(fnAllTableSumm)-1,ncol=ncol(fnAllTableSumm)-1)
    fnMatrix<-sweep(fnAllTableSumm[c(seq(nrow(fnAllTableSumm)-1)),c(seq(ncol(fnAllTableSumm)-1))],MARGIN=1,fnAllTableSumm[c(seq(nrow(fnAllTableSumm)-1)),ncol(fnAllTableSumm)],'*')
    fnMatrix<-sweep(fnMatrix,MARGIN=2,fnColSumm,'/')
    fnColSumm<-c(fnColSumm,sum(fnColSumm))
    fnAllTableSumm["1",]<-fnColSumm
    write.table(fnAllTableSumm,file=fnFileName,sep="\t",quote=FALSE,row.names=F)
    return(fnMatrix)
}

VennDiag<-function(fnSetList,fnOutputPath,fnOutputFileName)
{
    print("*************************  Running Venn Diagram  *******************")
    fnMethodToPrint<-paste("VennDiag(fnSetList,",fnOutputPath,",",fnOutputFileName,")",collapse="",sep="")
    print(fnMethodToPrint)
    
    if(length(unlist(fnSetList)>0))
    {
        #fnColor<-c("yellow","red","green","seagreen3","purple","orchid3")
       fnColor<-c("yellow","red","#56B4E9","green","purple","orchid3")
       fnOutputFileNameVenn<-paste(fnOutputPath,"/",fnOutputFileName,collapse="",sep = "")
       if(is.null(names(fnSetList))){
           names(fnSetList)<-paste("Set",1:length(fnSetList),sep="")}
       fnDataForVenn<-Filter(length,fnSetList)
       fnDataForVenn<-fnDataForVenn[names(sort(sapply(fnDataForVenn,length)))]
       ven.plot<-venn.diagram(fnDataForVenn,fill = fnColor[1:length(fnDataForVenn)],alpha =rep(0.5,length(fnDataForVenn)), cex = 2,cat.fontface = 4,lty =1,filename=paste(fnOutputFileNameVenn,".tiff",collapse="",sep = ""),main="",sub=fnOutputFileName,euler.d=TRUE,scaled=TRUE,main.fontface=0,sub.fontface=1,sep.dist=0,offset=0.7)
       ven.plot<-venn.diagram(fnDataForVenn,fill = fnColor[1:length(fnDataForVenn)],
                              imagetype="png",
                              height=128,
                              width=128,
                              resolution=20,
                              alpha =rep(0.5,length(fnDataForVenn)), 
                              cex = 2,cat.fontface = 4,lty =1,
                              filename=paste(fnOutputFileNameVenn,".png",collapse="",sep = ""),
                              main="",
                              sub=fnOutputFileName,
                              euler.d=TRUE,
                              scaled=TRUE,
                              main.fontface=0,
                              sub.fontface=1,
                              sep.dist=0,offset=0.7)
       #ven.plot<-venn.diagram(fnDataForVenn,fill = fnColor[1:length(fnDataForVenn)],alpha =rep(0.5,length(fnDataForVenn)), cex = 2,cat.fontface = 4,lty =1,filename=paste(fnOutputFileNameVenn,".tiff",collapse="",sep = ""),main="",euler.d=TRUE,scaled=TRUE,main.fontface=0,sub.fontface=1,sep.dist=0,offset=0.7,cat.dist=0.05)
       fnUniverse<-sort(unique(unlist(fnDataForVenn)))
       fnTableInters<- matrix(0, nrow=length(fnUniverse), ncol=length(fnSetList))
       rownames(fnTableInters)<-fnUniverse
       colnames(fnTableInters) <- names(fnSetList)
       for (j in 1:length(fnUniverse)) {
           for(i in 1:length(fnSetList)){
               fnTableInters[j,i] <- fnUniverse[j] %in% fnSetList[[i]]}
       }
       fnTableInters<-fnTableInters[order(rowSums(fnTableInters),decreasing=TRUE),]
       fnMethodWeightMatrix<-IntersectionSummary(fnTableInters,paste(fnOutputFileNameVenn,"_IntersectSummary.txt",collapse="",sep = ""))
       #fnGeneWeight<-as.double(rowSums(sweep(fnTableInters,MARGIN=2,fnMethodWeight,'*'))/sum(fnMethodWeight))
       #print(fnGeneWeight)
       #fnTableInters$Weight<-fnGeneWeight
       write.table(fnMethodWeightMatrix,file=paste(fnOutputFileNameVenn,"_matrixWeight.txt",collapse="",sep = ""),row.names=F, sep="\t",quote=FALSE)
       write.table(fnTableInters, file=paste(fnOutputFileNameVenn,"_table.txt",collapse="",sep = ""), sep="\t",quote=FALSE)
       UpSetPlot(fnSetList,fnOutputPath,fnOutputFileName)
       return(1)
    }
}

UpSetPlot<-function(fnSetList,fnOutputPath,fnOutputFileName)
{
    print("*************************  Running UpSetPlot  *******************")
    fnMethodToPrint<-paste("UpSetPlot(fnSetList,",fnOutputPath,",",fnOutputFileName,")",collapse="",sep="")
    print(fnMethodToPrint)
    fnColor<-c("yellow","red","#56B4E9","green","purple","orchid3")
    
    if(length(unlist(fnSetList)>0))
    {
        fnOutputFileNameUpSet<-paste(fnOutputPath,"/",fnOutputFileName,"_upsetR",collapse="",sep = "")
        if(is.null(names(fnSetList))){
            names(fnSetList)<-paste("Set",1:length(fnSetList),sep="")}
        fnDE_gns<-Filter(length,fnSetList)
        fnDE_gns <- UpSetR::fromList(fnDE_gns)
        pdf(paste(fnOutputFileNameUpSet,".pdf",collapse="",sep=""),onefile=FALSE)
        png(paste(fnOutputFileNameUpSet,".png",collapse="",sep=""),height=128, width=128,res=20)
        for(fnDevice in dev.list())
        {
            dev.set(fnDevice)
            print(UpSetR::upset(fnDE_gns,order.by = "freq",sets.bar.color = rev(fnColor[1:length(fnSetList)])))
        }
        #ggsave(fnOutputFileNameUpSet,UpSetR::upset(fnDE_gns,order.by = "freq",sets.bar.color = "#56B4E9"))
        graphics.off()
    }

}
