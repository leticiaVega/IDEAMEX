# source("~/bin/DifferentialExp/HeatmapPlots.r")
heatmapDEGenes<-function(fnDEListGenes,fnCountTable,fnOutputPath,fnConditions,fnNameSufix)
{
    if(!("package:ComplexHeatmap" %in% search()))
    {
        library(ComplexHeatmap)
    }
    if(!("package:DESeq2" %in% search()))
    {
        library(DESeq2)
    }
    if(!("package:circlize" %in% search()))
    {
        library(circlize)
    }
    fnDEGenes<-scan(fnDEListGenes,what=as.character())
    fnNumberOfGenes<-length(fnDEGenes)
    if(fnNumberOfGenes>200)
    {
        fnDEGenes<-fnDEGenes[1:200]
        fnNumberOfGenes<-200
    }
    fnGenesFontSize<-7 - ((fnNumberOfGenes - 1) %/% 50)

    ####  Iicializacion de variables
    fnSamplesName=factor(sub("_[a-zA-Z0-9]+$","",colnames(fnCountTable)))
    fnConditionsNames<-paste(fnConditions[1],"vs",fnConditions[2],collapse="",sep = "")
    fnFileName<-paste(fnOutputPath,"/",fnConditionsNames,fnNameSufix,collapse="",sep = "")

    ####  Initializacion de un objecto de tipo DESeq
    fnColData<-data.frame(condition=fnSamplesName,row.names=colnames(fnCountTable))
    fnDds <- try(DESeqDataSetFromMatrix(countData = as.matrix(fnCountTable),colData = fnColData,design = ~ condition),silent=TRUE)
    if(!(is(fnDds,"try-error")))
    {
        fnPlotFileName<-paste(fnFileName,"_heatmap",collapse="",sep = "")
        pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
        png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
        fnRld = rlog(fnDds, blind = FALSE)
        fnGenes<-assay(fnRld)[fnDEGenes,]
        fnZ.mat <- t(scale(t(fnGenes), center=TRUE, scale=TRUE))
        # colour palette
        fnMyPalette <- c("red3", "ivory", "blue3")
        fnMyRamp = colorRamp2(c(-2, 0, 2), fnMyPalette)
        for(fnDevice in dev.list())
        {
            dev.set(fnDevice)
            print(Heatmap(fnZ.mat, name = "z-score", col = fnMyRamp, show_row_name = TRUE, row_names_gp = gpar(fontsize = fnGenesFontSize),column_title = fnConditionsNames))
            #Heatmap(fnZ.mat, name = "z-score", col = fnMyRamp, show_row_name = FALSE, cluster_columns = FALSE)
        }
        #### Cierre del modo de guardado de graficos
        graphics.off()
        printOKMessage("       Heatmap .......................... OK")
    }
    else{
        printErrorMessage("      Objeto Dds .......................... Failed")
        return(fnDds)
    }
}





