# source("~/bin/DifferentialExp/RunGetInfoFromFile.r")
# RunDataAbundanceGraph("~/Temp/Flosan",c("edgeR","NOISeq"),"BOTTOMvsTOP")


RunGetInfoFromFile<-function(fnOutputPath,fnDEMethods,fnFileName)
{
    print("*************************  Running Get Info from File  *******************")
    fnVennFileName<-paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_table.txt",collapse="",sep = "")
    fnVennTable<-read.table(fnVennFileName,header=T,row.names=1,sep="\t")
    fnTop<-row.names(fnVennTable[rowSums(fnVennTable)==length(fnDEMethods),])
    for(i in fnDEMethods)
    {
        fnFileMethod<-paste(fnOutputPath,"/",i,"_Results/",fnFileName,".txt",collapse="",sep = "")
        fnTableMethod<-read.table(fnFileMethod,header=T,row.names=1,sep="\t")
        write.table(fnTableMethod[fnTop,],paste(fnOutputPath,"/",i,"_Results/",fnFileName,"_Intersect.txt",collapse="",sep = ""),sep="\t",quote=F,col.names=T,row.names=T)
    }
    write(fnTop,file=paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_Intesrsect_TOP_IDs.txt",collapse="",sep = ""))
    write(row.names(fnVennTable),file=paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_Union_TOP_IDs.txt",collapse="",sep = ""))
}

RunDataAbundance<-function(fnOutputPath,fnDEMethods,fnFileName,fnList)
{
    print("*************************  Running Data Abundance  *******************")
    ####  Cargando la biblioteca
    if(!("package:corrplot" %in% search()))
    {
        library(corrplot)
    }
    fnNewNames<-c()
    fnAllAbundance<-data.frame(row.names=fnList)
    fnAllAbundance<-fnAllAbundance[order(row.names(fnAllAbundance)),]
    for(i in fnDEMethods)
    {
        fnFileMethod<-paste(fnOutputPath,"/",i,"_Results/",fnFileName,"_Abundances.txt",collapse="",sep = "")
        fnTableMethod<-read.table(fnFileMethod,header=T,row.names=1,sep="\t",stringsAsFactors =FALSE)
        fnAbssentnames<-setdiff(row.names(fnAllAbundance),row.names(fnTableMethod))
        if(length(fnAbssentnames)>0)
        {
            fnFixedAbbss<-data.frame(row.names = fnAbssentnames,matrix(data = 0,nrow = length(fnAbssentnames), ncol=ncol(fnTableMethod)))
            names(fnFixedAbbss)<-names(fnTableMethod)
            fnTableMethod<-rbind(fnTableMethod,fnFixedAbbss)
        }
        else{
            fnTableMethod<-fnTableMethod[fnList,]
        }
        fnTableMethod<-fnTableMethod[order(row.names(fnTableMethod)),]
        fnAllAbundance<-cbind(fnAllAbundance,fnTableMethod)
        fnNewNames<-c(fnNewNames,paste(i,names(fnTableMethod),sep="_"))
    }
    fnAllCondsWithoutindex=factor(sub("_[a-zA-Z0-9]+$","",colnames(fnAllAbundance)))
    names(fnAllAbundance)<-fnNewNames
    fnCorrPlotName<-paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_AbundanceCorrelation",collapse="",sep = "")
    pdf(paste(fnCorrPlotName,".pdf",collapse="",sep=""))
    for(i in levels(fnAllCondsWithoutindex))
    {
        fnSubsetCondtion<-fnAllAbundance[,grep(i, names(fnAllAbundance), value=TRUE)]
        fnCorrelation = cor(fnSubsetCondtion,method="spearman")
        corrplot(fnCorrelation, method = "number",type = "upper",number.digits=4,tl.cex = 0.6, tl.srt = 45,number.cex=0.6, mar = c(2, 2, 2, 2),xpd = TRUE,main=i)
    }
    png(paste(fnCorrPlotName,".png",collapse="",sep=""),height=128, width=128,res=20)
    corrplot(fnCorrelation, method = "number",type = "upper",number.digits=4,tl.cex = 0.6, tl.srt = 45,number.cex=0.6, mar = c(2, 2, 2, 2),xpd = TRUE,main=levels(fnAllCondsWithoutindex)[0])
    graphics.off()
    write.table(fnAllAbundance,paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_AbundanceTable.txt",collapse="",sep = ""),sep="\t",quote=F,col.names=T,row.names=T)
}

RunDataLogFC<-function(fnOutputPath,fnDEMethods,fnFileName,fnList)
{
    print("*************************  Running Data LogFC  *******************")
    ####  Cargando la biblioteca
    if(!("package:corrplot" %in% search()))
    {
        library(corrplot)
    }
    fnNewNames<-c()
    fnAllLogFC<-data.frame(row.names=fnList)
    fnAllLogFC<-fnAllLogFC[order(row.names(fnAllLogFC)),]
    for(i in fnDEMethods)
    {
        fnFileMethod<-paste(fnOutputPath,"/",i,"_Results/",fnFileName,"_logFC.txt",collapse="",sep = "")
        fnTableMethod<-read.table(fnFileMethod,header=T,row.names=1,sep="\t",stringsAsFactors =FALSE)
        fnAbssentnames<-setdiff(row.names(fnAllLogFC),row.names(fnTableMethod))
        if(length(fnAbssentnames)>0)
        {
            fnFixedAbbss<-data.frame(row.names = fnAbssentnames,matrix(data = 0,nrow = length(fnAbssentnames), ncol=ncol(fnTableMethod)))
            names(fnFixedAbbss)<-names(fnTableMethod)
            fnTableMethod<-rbind(fnTableMethod,fnFixedAbbss)
        }
        else{
            fnTableMethod<-fnTableMethod[fnList,,drop=FALSE]
        }
        fnTableMethod<-fnTableMethod[order(row.names(fnTableMethod)),]
        fnAllLogFC<-cbind(fnAllLogFC,fnTableMethod)
    }
    names(fnAllLogFC)<-fnDEMethods
    fnAllLogFC<-na.omit(fnAllLogFC)
    fnCorrPlotName<-paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_logFCCorrelation",collapse="",sep = "")
    fnCorrelation = cor(fnAllLogFC,method="spearman")
    pdf(paste(fnCorrPlotName,".pdf",collapse="",sep=""))
    corrplot(fnCorrelation, method = "number",type = "upper",number.digits=4,tl.cex = 0.6, tl.srt = 45,number.cex=0.6, mar = c(2, 2, 2, 2),xpd = TRUE,main=paste("logFC correlation",fnFileName))
    png(paste(fnCorrPlotName,".png",collapse="",sep=""),height=128, width=128,res=20)
    corrplot(fnCorrelation, method = "number",type = "upper",number.digits=4,tl.cex = 0.6, tl.srt = 45,number.cex=0.6, mar = c(2, 2, 2, 2),xpd = TRUE,main=paste("logFC correlation",fnFileName))
    graphics.off()
    write.table(fnAllLogFC,paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_logFCTable.txt",collapse="",sep = ""),sep="\t",quote=F,col.names=T,row.names=T)
}

RunDataPval<-function(fnOutputPath,fnDEMethods,fnFileName,fnList)
{
    print("*************************  Running Data Pval  *******************")
    ####  Cargando la biblioteca
    if(!("package:corrplot" %in% search()))
    {
        library(corrplot)
    }
    fnNewNames<-c()
    fnAllPval<-data.frame(row.names=fnList)
    fnAllPval<-fnAllPval[order(row.names(fnAllPval)),]
    for(i in fnDEMethods)
    {
        fnFileMethod<-paste(fnOutputPath,"/",i,"_Results/",fnFileName,"_pval.txt",collapse="",sep = "")
        fnTableMethod<-read.table(fnFileMethod,header=T,row.names=1,sep="\t",stringsAsFactors =FALSE)
        fnAbssentnames<-setdiff(row.names(fnAllPval),row.names(fnTableMethod))
        if(length(fnAbssentnames)>0)
        {
            fnFixedAbbss<-data.frame(row.names = fnAbssentnames,matrix(data = 0,nrow = length(fnAbssentnames), ncol=ncol(fnTableMethod)))
            names(fnFixedAbbss)<-names(fnTableMethod)
            fnTableMethod<-rbind(fnTableMethod,fnFixedAbbss)
        }
        else{
            fnTableMethod<-fnTableMethod[fnList,,drop=FALSE]
        }
        fnTableMethod<-fnTableMethod[order(row.names(fnTableMethod)),]
        fnAllPval<-cbind(fnAllPval,fnTableMethod)
    }
    names(fnAllPval)<-fnDEMethods
    fnAllPval<-na.omit(fnAllPval)
    fnCorrPlotName<-paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_PvalCorrelation",collapse="",sep = "")
    fnCorrelation = cor(fnAllPval,method="spearman")
    pdf(paste(fnCorrPlotName,".pdf",collapse="",sep=""))
    corrplot(fnCorrelation, method = "number",type = "upper",number.digits=4,tl.cex = 0.6, tl.srt = 45,number.cex=0.6, mar = c(2, 2, 2, 2),xpd = TRUE,main=paste("Pval correlation",fnFileName))
    png(paste(fnCorrPlotName,".png",collapse="",sep=""),height=128, width=128,res=20)
    corrplot(fnCorrelation, method = "number",type = "upper",number.digits=4,tl.cex = 0.6, tl.srt = 45,number.cex=0.6, mar = c(2, 2, 2, 2),xpd = TRUE,main=paste("Pval correlation",fnFileName))
    graphics.off()
    write.table(fnAllPval,paste(fnOutputPath,"/VennDiagram_Results/",fnFileName,"_PvalTable.txt",collapse="",sep = ""),sep="\t",quote=F,col.names=T,row.names=T)
}

