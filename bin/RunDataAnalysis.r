# source(~/bin/DifferentialExp/RunDataAnalysis.r)
# source(~/bin/DifferentialExp/RunPrintMessage.r)

dataframe2stack<-function(fnTable)
{
    fnNewTable<-stack(fnTable,select=colnames(fnTable))
    fnNewTable$Condition<-sub("_[a-zA-Z0-9]+$","",fnNewTable$ind)
    fnNewTable<-fnNewTable[,c(2,1,3)]
    names(fnNewTable)<-c("Samples","value","Condition")
    return(fnNewTable)
}

callMDS<-function(fnTable,fnPlotFileName,fnCondition,fnNormalization=FALSE)
{
    pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
    png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
    fnDge<-DGEList(counts=fnTable, group=factor(fnCondition))
    if(fnNormalization){fnDge=calcNormFactors(fnDge)}
    fnColors=as.numeric(fnDge$samples$group)+1
    for(fnDevice in dev.list())
    {
        dev.set(fnDevice)
        plotMDS(fnDge,col=fnColors)
    }
    graphics.off()
}

callPCA<-function(fnTable,fnPlotFileName,fnCondition,fnNormalization=FALSE)
{
    fnCds<-newCountDataSet(fnTable,fnCondition)
    if(fnNormalization){fnCds <- DESeq::estimateSizeFactors(fnCds)}
    else{fnCds$sizeFactor<-rep(1,length(fnCondition))}
    fnCds <- DESeq::estimateSizeFactors(fnCds)
    if(length(fnCondition)==length(levels(factor(fnCondition)))){fnCdsB<-DESeq::estimateDispersions(fnCds,method="blind",fitType="local",sharingMode="fit-only")}
    else{fnCdsB<-DESeq::estimateDispersions(fnCds,method="pooled",fitType="local",sharingMode="fit-only")} #,sharingMode="maximum"
    fnVsd<-DESeq::varianceStabilizingTransformation(fnCdsB)
    pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
    png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
    for(fnDevice in dev.list())
    {
        dev.set(fnDevice)
        print(DESeq::plotPCA(fnVsd,intgroup="condition"))
    }
    graphics.off()
}

plotCPM<-function(fnTable,fnPlotFileName)
{
    fnCPMsTable = cpm(fnTable)
    fnCPMHistogram<-data.frame()
    fnCPMHistogram<- rbind(fnCPMHistogram,colSums(fnCPMsTable>10))
    fnCPMHistogram<- rbind(fnCPMHistogram,colSums(fnCPMsTable>5 & fnCPMsTable<= 10))
    fnCPMHistogram<- rbind(fnCPMHistogram,colSums(fnCPMsTable>=3 & fnCPMsTable<= 5))
    fnCPMHistogram<- rbind(fnCPMHistogram,colSums(fnCPMsTable>=2 & fnCPMsTable<3))
    fnCPMHistogram<- rbind(fnCPMHistogram,colSums(fnCPMsTable>=1 & fnCPMsTable<2))
    fnCPMHistogram<- rbind(fnCPMHistogram,colSums(fnCPMsTable < 1))
    names(fnCPMHistogram)<-colnames(fnCPMsTable)
    rownames(fnCPMHistogram)<-c("CPM>10","5< CPM <= 10","2 < CPM <= 5","CPM=2","CPM=1","CPM=0")
    fnCPMHistogram<-as.matrix(fnCPMHistogram,byrow=T)
    pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
    png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
    for(fnDevice in dev.list())
    {
        dev.set(fnDevice)
        par(mar=c(10, 4.1, 4.1, 7), xpd=TRUE)
        barplot(fnCPMHistogram, col=c(8,colors()[133],colors()[128],colors()[90],6,colors()[657]),las=3,cex.axis = 1.0,cex.names=0.7,width=2)
        legend("topright",inset=c(-0.25,0), fill=c(8,colors()[133],colors()[128],colors()[90],6,colors()[657]), legend=rownames(fnCPMHistogram),cex=0.7)
    }
    graphics.off()
}

callBoxPlot<-function(fnTable,fnPlotFileName,fnCondition,fnYlab)
{
    fnColors=as.numeric(factor(fnCondition))+1
    pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
    png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
    for(fnDevice in dev.list())
    {
        dev.set(fnDevice)
        boxplot(fnTable, col=fnColors,las=3, ylab=fnYlab,par(mar = c(8, 5, 4, 2)+ 0.1),outline=FALSE) # ,cex.axis = 1.0
    }
    graphics.off()
}

callDensityPlot<-function(fnTable,fnPlotFileName)
{
    fnDf = dataframe2stack(fnTable)
    fnDf = data.frame(fnDf, Condition = sub("_[a-zA-Z0-9]+$","",fnDf$Samples))
    pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
    png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
    for(fnDevice in dev.list())
    {
        dev.set(fnDevice)
        print(ggplot(fnDf, aes(x=value,colour = Samples, fill = Samples)) + ylim(c(0, 0.25)) + geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) + xlab(expression(log[2](count + 1))))
    }
    graphics.off()
}

dataAnalysis<-function(fnCountTable,fnCondition,fnOutputPath)
{
   library(ggplot2)
   library(DESeq)
   library(DESeq2)
   library(edgeR)
   library(NOISeq)
   library(RColorBrewer)

   fnPseudoCountTable<-log2(fnCountTable + 1)
   fnExpDesign<-data.frame(fnSamples = colnames(fnCountTable),fnFactor =fnCondition)
   fnMyData<-readData( data=fnCountTable, factors=fnExpDesign)
   fnNormalizedCountTable<-data.frame(tmm(assayData(fnMyData)$exprs, k=0.5, lc = 0))
   
   
   fnMethods<-c("callBoxPlot(fnPseudoCountTable,fnPlotFileName,fnCondition,\"log[2](count + 1)\")","callDensityPlot(fnPseudoCountTable,fnPlotFileName)","plotCPM(fnCountTable,fnPlotFileName)","callPCA(fnCountTable,fnPlotFileName,fnCondition)","callMDS(fnCountTable,fnPlotFileName,fnCondition)","callBoxPlot(fnNormalizedCountTable,fnPlotFileName,fnCondition,\"TMM normalization\")","callPCA(fnCountTable,fnPlotFileName,fnCondition,fnNormalization=TRUE)","callMDS(fnCountTable,fnPlotFileName,fnCondition,fnNormalization=TRUE)")
   fnExt<-c("BoxPlotLog2","DensitiesPlot","CPMPlot","PCAPlot","MDSPlot","BoxPlotNorm","PCAPlotNorm","MDSPlotNorm")
   
   for(i in 1:length(fnMethods))
   {
       fnPlotFileName<-paste(fnOutputPath,"_",fnExt[i],collapse="",sep = "")
       if(is(try(eval(parse(text=fnMethods[i])),silent=TRUE),"try-error")){
           printErrorMessage(paste("    ",fnExt[i],"    .......................... Failed"))
           if(!(is.null(dev.list()))){
               graphics.off()}
           fnFilestoRemove<-list.files(path=dirname(fnPlotFileName),pattern=basename(fnPlotFileName),full.names=TRUE)
           if(length(fnFilestoRemove)>0){
               file.remove(fnFilestoRemove)}
       }
       else{
           printOKMessage(paste("   ",fnExt[i],"    .......................... OK"))}
   }
}

#plotAbundance(fnAbundance,fnAbundanceNormalized,fnConditions,fnMain,fnOutputPath)
#{
#    fnAllData<-cbind(fnAbundance,fnAbundanceNormalized)
#    par(mfrow=c(length(fnConditions),ncols(fnAbunda)/2))
#    for(i in seq(ncols(fnAbundance)))
#    {
#        j<-i+ncols(fnAbundance)
#        if(max(fnAllData[,i])>max(fnAllData[,j]))
#        {fnSortCol<-i}
#        else{fnSortCol<-j}
#        plot(fnAbundance[,i],col='blue',main="Frecuencia por ciclos",xlab="Ciclos",type="l",ylim=c(0,1),ylab="A+T+C+G",yaxt="n")
#        axis(2,at=c(0,0.25,0.5,0.75,1),las=2)
#        lines(fnAbundanceNormalized[,i],col='red')
#
#
#
#
#        plot(fnAbundance[,1],
#    }
#}






