RunDESeq2<- function(fnCountTable,fnOutputPath,TOP=FALSE,fnUmbral=0.01,fnUmbralFoldChange=1,fnMAPlot=TRUE,fnBatch=c(),fnConditions)
{
   print("*************************  Running DESeq2  *************************")
   fnMethodToPrint<-paste("RunDESeq2(fnCounTable,",fnOutputPath,",TOP=",TOP,",fnUmbral=",fnUmbral,",fnUmbralFoldChange=",fnUmbralFoldChange,",fnMAPlot=",fnMAPlot,",fnConditions=c(",fnConditions[1],",",fnConditions[2],")",")",collapse="",sep="")
   print(fnMethodToPrint)
   
   ####  Cargando la biblioteca
   if(!("package:DESeq2" %in% search()))
   {
       library(DESeq2)
   }
   ####  Iicializacion de variables
   fnSamplesName=factor(sub("_[a-zA-Z0-9]+$","",colnames(fnCountTable)))
   fnConditionsNames<-paste(fnConditions[1],"vs",fnConditions[2],collapse="",sep = "")
   fnFileName<-paste(fnOutputPath,"/",fnConditionsNames,collapse="",sep = "")
   fnOutputFileNameTop<-paste(fnFileName,"_TOP.txt",collapse="",sep = "")
   fnOutputFileName<-paste(fnFileName,".txt",collapse="",sep = "")
   print("############")
   print(paste("Samples: ",fnConditionsNames))
   print("############")
   
   ####  Initializacion de un objecto de tipo DESeq
   fnColData<-data.frame(condition=fnSamplesName,row.names=colnames(fnCountTable))
   if(length(fnBatch))
   {
       #printOKMessage("      Objeto Dds .......................... OK")
       printOKMessage("      Batch effects .......................... OK")
       fnColData$fnBatch<-factor(fnBatch)
       fnDds <- try(DESeqDataSetFromMatrix(countData = as.matrix(fnCountTable),colData = fnColData,design = ~ fnBatch + condition),silent=TRUE)
   }
   else{
       fnDds <- try(DESeqDataSetFromMatrix(countData = as.matrix(fnCountTable),colData = fnColData,design = ~ condition),silent=TRUE)
   }
   if(!(is(fnDds,"try-error")))
   {
       printOKMessage("      Objeto Dds .......................... OK")
       fnDds$condition <- factor(fnDds$condition, levels=c(fnConditions[1],fnConditions[2]))
       ####  Grafica de agrupamiento de los datos
       if(length(fnColData$condition)>2)
       {
           fnPlotFileName<-paste(fnFileName,"_plotPCA",collapse="",sep = "")
           pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
           png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
           fnRld<-rlogTransformation(fnDds, blind=FALSE)
           for(fnDevice in dev.list())
           {
               dev.set(fnDevice)
               print(DESeq2::plotPCA(fnRld, intgroup="condition"))
           }
           #### Cierre del modo de guardado de graficos
           graphics.off()
           printOKMessage("      PCA plot .......................... OK")
       }
       #### Calculo de la Expresion diferencial
       fnDds<-DESeq2::DESeq(fnDds)
       fnRes<-DESeq2::results(fnDds)
       printOKMessage("      Differential expression estimation.......................... OK")
       
       ####  Ordenando datos por nombre de genes para incorporar informacion adicional
       fnRes <- fnRes[order(row.names(fnRes)),]

       ####  Calculo de la abundancia normalizada
       fnAbundance<-cbind(DESeq2::counts(fnDds),DESeq2::counts(fnDds,normalized=TRUE))
       fnAbundance<-fnAbundance[order(row.names(fnAbundance)),]
       fnAbundanceFileName<-paste(fnFileName,"_Abundances.txt",collapse="",sep = "")
       write.table(DESeq2::counts(fnDds,normalized=TRUE), file=fnAbundanceFileName, sep="\t",quote=FALSE)
       printOKMessage("      Abundance estimation .......................... OK")
       
       fnFinalTab<-cbind(fnRes[,],fnAbundance[,])
       fnFinalTab <- fnFinalTab[order(fnFinalTab$padj,decreasing=FALSE),]
       fnPvalFileName<-paste(fnFileName,"_pval.txt",collapse="",sep = "")
       write.table(fnFinalTab[,"padj",drop=FALSE], file=fnPvalFileName, sep="\t",quote=FALSE)
       fnlogFCFileName<-paste(fnFileName,"_logFC.txt",collapse="",sep = "")
       write.table(fnFinalTab[,"log2FoldChange",drop=FALSE], file=fnlogFCFileName, sep="\t",quote=FALSE)

       ####  Calculo de la certeza y la regulacion
       fnFinalTab$Certeza<-round((1-fnFinalTab$padj)*100,2)
       fnFinalTab$DE <- ifelse(fnFinalTab$log2FoldChange <0, paste("UP_",fnConditions[1],"_DOWN_",fnConditions[2],collapse="",sep = ""),paste("DOWN_",fnConditions[1],"_UP_",fnConditions[2],collapse="",sep = ""))
       printOKMessage("      Results table .......................... OK")
    
       ####  Generando la grafica de MA
       fnPlotFileName<-paste(fnFileName,"_plotMA",collapse="",sep = "")
       pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
       png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
       for(fnDevice in dev.list())
       {
           dev.set(fnDevice)
           DESeq::plotMA(data.frame(fnRes), ylim=c(-4, 4), main=fnConditionsNames,col = ifelse((fnRes$padj<=fnUmbral) & (abs(fnRes$log2FoldChange)>=fnUmbralFoldChange), "red","black"),linecol = "blue",xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change),log = "x",cex=0.45)
           abline(h =c(-fnUmbralFoldChange, fnUmbralFoldChange), col = "blue") #c(-2, 2)
       }
       graphics.off()
       printOKMessage("      Results table .......................... OK")
       
       ####  Guardado de los datos en arhivo
       if(TOP)
       {
           fnTopRows<-which((fnFinalTab$padj<=fnUmbral) & (abs(fnFinalTab$log2FoldChange)>=fnUmbralFoldChange))
           write.table(fnFinalTab[fnTopRows,],file=fnOutputFileNameTop, sep="\t",quote=FALSE,row.names=TRUE)
           fnTopName<-row.names(fnFinalTab[fnTopRows,])
           printOKMessage("      Top file .......................... OK")
       }
       else
       {
           fnTopName=c()
       }
       write.table(fnFinalTab, file=fnOutputFileName, sep="\t",quote=FALSE,row.names=TRUE)
       printOKMessage("      Final Results table to file .......................... OK")
       return(fnTopName)
   }
   else{
       printErrorMessage("      Objeto Dds .......................... Failed")
       return(fnDds)
   }
   
}
