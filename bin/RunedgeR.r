# source("~/bin/DifferentialExp/RunedgeR.r")
# source("~/bin/DifferentialExp/RunPrintMessage.r")

evalDispTestWithOutRep<-function(fnDge,fnConditions)
{
    fnDispersion = 0.4
    fnEt = exactTest(fnDge, dispersion=fnDispersion , pair=c(fnConditions[1],fnConditions[2]))
    printOKMessage("      Differential expression estimation.......................... OK")
    return(fnEt)
}

evalDispTestWithRep<-function(fnDge,fnConditions)
{
    ####  Calculo de la dispersion de los datos
    fnDge = estimateDisp(fnDge)
    printOKMessage("      Dispersion estimation .......................... OK")
    ####  Calculo de la Expresion diferencial
    fnEt = exactTest(fnDge,pair=c(fnConditions[1],fnConditions[2]))
    printOKMessage("      Differential expression estimation.......................... OK")
    return(fnEt)
}

evalDispTestWithBatch<-function(fnDge,fnConditionsRef,fnSamplesName,fnBatch)
{
    ####  Calculo de la dispersion de los datos
    fnSamplesName<-relevel(fnSamplesName,ref=fnConditionsRef)
    fnDesign<-model.matrix(~fnBatch+fnSamplesName)
    fnDge = estimateDisp(fnDge, fnDesign, robust = TRUE)
    fnFit <- glmQLFit(fnDge, fnDesign, robust = TRUE)
    fnQlf <- glmQLFTest(fnFit)
    printOKMessage("      Differential expression estimation.......................... OK")
    return(fnQlf)
}


RunedgeR<-function(fnCountTable,fnOutputPath,TOP=FALSE,fnUmbral=0.01,fnUmbralFoldChange=1,fnSmearPlot=TRUE,fnVolcanoPlot=TRUE,fnBatch=c(),fnConditions)
{
   print("*************************  Running edgeR  *************************")
   fnMethodToPrint<-paste("RunedgeR(fnCounTable,",fnOutputPath,",TOP=",TOP,",fnUmbral=",fnUmbral,",fnUmbralFoldChange=",fnUmbralFoldChange,",fnSmearPlot=",fnSmearPlot,",fnVolcanoPlot=",fnVolcanoPlot,",fnConditions=c(",fnConditions[1],",",fnConditions[2],")",")",collapse="",sep="")
   print(fnMethodToPrint)
   ####  Cargando la biblioteca
   if(!("package:edgeR" %in% search()))
   {
       library(edgeR)
   }
   ####  Iicializacion de variables
   fnSamplesName=factor(sub("_[a-zA-Z0-9]+$","",colnames(fnCountTable)))
   fnConditionsNames<-paste(fnConditions[1],"vs",fnConditions[2],collapse="",sep = "")
   print("############")
   print(paste("Samples: ",fnConditionsNames))
   print("############")
   fnFileName<-paste(fnOutputPath,"/",fnConditionsNames,collapse="",sep = "")
   fnPlotFileName<-paste(fnFileName,"_plotMDS",collapse="",sep = "")
   fnOutputFileName<-paste(fnFileName,".txt",collapse="",sep = "")
   fnOutputFileNameTop<-paste(fnFileName,"_TOP.txt",collapse="",sep = "")
   
   ### Composición del objeto DGEList
   fnDge<-try(DGEList(counts=fnCountTable, group=fnSamplesName),silent=TRUE)
   if(!(is(fnDge,"try-error")))
   {
       printOKMessage("      Objeto DGEList .......................... OK")
       ####  Normalizacion de los datos
       fnDge=calcNormFactors(fnDge)
       printOKMessage("      Normalizacion .......................... OK")
       ####  Grafica de agrupamiento de los datos
       if(length(fnDge$samples$group)>2)
       {
           fnColors=as.numeric(fnDge$samples$group)
           pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
           png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
           for(fnDevice in dev.list())
           {
               dev.set(fnDevice)
               plotMDS(fnDge,col=fnColors)
           }
           #### Cierre del modo de guardado de graficos
           graphics.off()
           printOKMessage("      MDS Plot .......................... OK")
       }
       if(length(fnDge$samples$group)<=2)
       {
           fnEt<-evalDispTestWithOutRep(fnDge,fnConditions)
       }
       else if(length(fnBatch))
       {
           fnEt<-evalDispTestWithBatch(fnDge,fnConditions[1],fnSamplesName,factor(fnBatch))
       }
       else
       {
           fnEt<-evalDispTestWithRep(fnDge,fnConditions)
       }
       
       ####  Obtencion de la tabla de resultados
       fnDeTab=topTags(fnEt,n=Inf)$table
       
       ####  Ordenando datos por nombre de genes para incorporar informacion adicional
       fnDeTab <- fnDeTab[order(row.names(fnDeTab)),]
       
       ####  Calculo de la abundancia normalizada
       fnAbundance<-cbind(fnDge$counts[,],cpm(fnDge,normalized.lib.size=T))
       fnAbundance<-fnAbundance[order(row.names(fnAbundance)),]
       fnAbundanceFileName<-paste(fnFileName,"_Abundances.txt",collapse="",sep = "")
       write.table(cpm(fnDge,normalized.lib.size=T), file=fnAbundanceFileName, sep="\t",quote=FALSE)
       printOKMessage("      Abundance estimation .......................... OK")
       
       fnFinalTab<-cbind(fnDeTab[,],fnAbundance[,])
       fnFinalTab <- fnFinalTab[order(fnFinalTab$FDR),]
       fnPvalFileName<-paste(fnFileName,"_pval.txt",collapse="",sep = "")
       write.table(fnFinalTab[,"FDR",drop=FALSE], file=fnPvalFileName, sep="\t",quote=FALSE)
       fnlogFCFileName<-paste(fnFileName,"_logFC.txt",collapse="",sep = "")
       write.table(fnFinalTab[,"logFC",drop=FALSE], file=fnlogFCFileName, sep="\t",quote=FALSE)
       
       ####  Calculo de la certeza y la regulacion
       fnFinalTab$DE<-round((1-fnFinalTab$FDR)*100,2)
       fnFinalTab$Regulation <- ifelse(fnFinalTab$logFC <0,paste("UP_",fnConditions[1],"_DOWN_",fnConditions[2],collapse="",sep = ""),paste("DOWN_",fnConditions[1],"_UP_",fnConditions[2],collapse="",sep = ""))
       printOKMessage("      Results table .......................... OK")
       
       ####  Generando la grafica de Smear
       if(fnSmearPlot)
       {
           deGenes = rownames(fnDeTab)[fnDeTab$FDR <= fnUmbral & abs(fnDeTab$logFC) >= fnUmbralFoldChange]
           fnPlotFileName<-paste(fnFileName,"_plotSmear",collapse="",sep = "")
           pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
           png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
           for(fnDevice in dev.list())
           {
               dev.set(fnDevice)
               plotSmear(fnDge, de.tags=deGenes,cex=0.5,main=fnConditionsNames,pair=c(fnConditions[1],fnConditions[2]),xlab="Log Concentration" , ylab="Log Fold-Change")
               abline(h =c(-fnUmbralFoldChange, fnUmbralFoldChange), col = "blue") #c(-2, 2)
           }
           graphics.off()
           printOKMessage("      Smear plot .......................... OK")
       }
       
       ####  Generando la grafica de Volcano
       if(fnVolcanoPlot)
       {
           fnPlotFileName<-paste(fnFileName,"_plotVolcano",collapse="",sep = "")
           pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
           png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
           fnLogF_log10 <- data.frame(logFC = fnFinalTab[, "logFC"], negLogPval = -log10(fnFinalTab[, "PValue"]))
           signGenes = (abs(fnLogF_log10$logFC) > fnUmbralFoldChange & fnLogF_log10$negLogPval > -log10(fnUmbral))
           for(fnDevice in dev.list())
           {
               dev.set(fnDevice)
               par(mar = c(5, 4, 4, 4))
               plot(fnLogF_log10, pch = 16, cex = 0.6,main=fnConditionsNames, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
               points(fnLogF_log10[signGenes, ], pch = 16, cex = 0.8, col = "red")
               abline(h = -log10(fnUmbral), col = "green3", lty = 2)
               abline(v = c(-fnUmbralFoldChange, fnUmbralFoldChange), col = "blue", lty = 2)
               mtext(paste("pval =", fnUmbral), side = 4, at = -log10(fnUmbral), cex = 0.8, line = 0.5, las = 1)
               mtext(c(paste("-", fnUmbralFoldChange, "fold"), paste("+", fnUmbralFoldChange, "fold")), side = 3, at = c(-fnUmbralFoldChange, fnUmbralFoldChange),cex = 0.8, line = 0.5)
           }
           graphics.off()
           printOKMessage("      Volcano plot .......................... OK")
       }
       ####  Guardado de los datos en archivo
       if(TOP)
       {
           fnTopRows<-which((fnFinalTab$FDR<=fnUmbral) & (abs(fnFinalTab$logFC) >= fnUmbralFoldChange))
           write.table(fnFinalTab[fnTopRows,],file=fnOutputFileNameTop, sep="\t",quote=FALSE)
           fnTopName<-row.names(fnFinalTab[fnTopRows,])
           printOKMessage("      Top file .......................... OK")
       }
       else
       {
           fnTopName=c()
       }
       write.table(fnFinalTab, file=fnOutputFileName, sep="\t",quote=FALSE)
       printOKMessage("      Final Results table to file .......................... OK")
       return(fnTopName)
   }
   else{
       printErrorMessage("      Objeto DGEList .......................... Failed")
       return(fnDge)
   }
}
