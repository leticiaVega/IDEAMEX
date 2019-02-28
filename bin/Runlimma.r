Runlimma<- function(fnCountTable,fnOutputPath,TOP=FALSE,fnUmbral=0.01,fnUmbralFoldChange=1,fnMDPlot=TRUE,fnBatch=c(),fnConditions)
{
   print("*************************  Running limma-Voom  *************************")
   fnMethodToPrint<-paste("Runlimma(fnCounTable,",fnOutputPath,",TOP=",TOP,",fnUmbral=",fnUmbral,",fnUmbralFoldChange=",fnUmbralFoldChange,",fnMDPlot=",fnMDPlot,",fnConditions=c(",fnConditions[1],",",fnConditions[2],")",")",collapse="",sep="")
   print(fnMethodToPrint)
   
   if(ncol(fnCountTable) > 2)
   {
       ####  Cargando la biblioteca
       if(!("package:limma" %in% search()))
       {
           library(limma)
       }
       ####  Inicializacion de variables
       fnSamplesName=factor(sub("_[a-zA-Z0-9]+$","",colnames(fnCountTable)),levels=c(fnConditions[1],fnConditions[2])) #
       fnConditionsNames<-paste(fnConditions[1],"vs",fnConditions[2],collapse="",sep = "")
       print("############")
       print(paste("Samples: ",fnConditionsNames))
       print("############")
       fnFileName<-paste(fnOutputPath,"/",fnConditionsNames,collapse="",sep = "")
       fnPlotFileName<-paste(fnFileName,"_plotMDS",collapse="",sep = "")
       fnOutputFileName<-paste(fnFileName,".txt",collapse="",sep = "")
       fnOutputFileNameTop<-paste(fnFileName,"_TOP.txt",collapse="",sep = "")
       if(length(fnBatch))
       {
           fnBatch<-factor(fnBatch)
           fnDesign<-model.matrix(~fnBatch+fnSamplesName)
           print("      Batch effect .......................... OK")
       }
       else{
           fnDesign <- model.matrix (~ fnSamplesName)
       }
       rownames(fnDesign) <- colnames(fnCountTable)
 
       ### Composición del objeto DGEList
       fnDge<-try(DGEList(counts=fnCountTable, group=fnSamplesName),silent=TRUE)
       if(!(is(fnDge,"try-error")))
       {
           printOKMessage("      Objeto DGEList .......................... OK")

           ####  Normalizacion de los datos
           fnDge=calcNormFactors(fnDge)
           printOKMessage("      Normalizacion .......................... OK")
           ####  Grafica de agrupamiento de los datos
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
           #### transform the count data to log2 -counts -per - million and estimate
           #### the mean - variance relationship , which is used to compute weights
           #### for each count -- this is supposed to make the read counts
           #### amenable to be used with linear models
           fnVoomTransformed <- voom(fnDge, fnDesign, plot = FALSE)
           printOKMessage("      Voom transform .......................... OK")
           ####  Ajustar un modelo lineal para cada gen.
           fnVoomedFitted <- lmFit ( fnVoomTransformed , design = fnDesign )
           printOKMessage("      lmfit estimation .......................... OK")
           ####  Calcular estadísticas t moderadas, estadísticas F moderadas y registro - probabilidades de expresión diferencial
           fnVoomedFitted <- eBayes(fnVoomedFitted)
           printOKMessage("      Differential expression estimation.......................... OK")
   
           ####  Obtencion de la tabla de resultados
           fnDeTab=topTable(fnVoomedFitted, coef = paste("fnSamplesName",fnConditions[2],sep="",collapse=""), number = Inf , adjust.method = "BH")
           ####  Ordenando datos por nombre de genes para incorporar informacion adicional
           fnDeTab <- fnDeTab[order(row.names(fnDeTab)),]
       
           ####  Calculo de la abundancia normalizada
           #fnAbundance<-cbind(fnDge$counts[,],cpm(fnDge, log=TRUE, prior.count=fnFilterValue)) Así estaba
           fnAbundance<-cbind(fnDge$counts[,],cpm(fnDge,normalized.lib.size=T))
           fnAbundance<-fnAbundance[order(row.names(fnAbundance)),]
           fnAbundanceFileName<-paste(fnFileName,"_Abundances.txt",collapse="",sep = "")
           write.table(cpm(fnDge,normalized.lib.size=T), file=fnAbundanceFileName, sep="\t",quote=FALSE)
           printOKMessage("      Abundance estimation .......................... OK")
       
           fnFinalTab<-cbind(fnDeTab[,],fnAbundance[,])
           fnFinalTab <- fnFinalTab[order(fnFinalTab$adj.P.Val),]
           fnPvalFileName<-paste(fnFileName,"_pval.txt",collapse="",sep = "")
           write.table(fnFinalTab[,"adj.P.Val",drop=FALSE], file=fnPvalFileName, sep="\t",quote=FALSE)
           fnlogFCFileName<-paste(fnFileName,"_logFC.txt",collapse="",sep = "")
           write.table(fnFinalTab[,"logFC",drop=FALSE], file=fnlogFCFileName, sep="\t",quote=FALSE)
       
           ####  Calculo de la certeza y la regulacion
           fnFinalTab$DE<-round((1-fnFinalTab$adj.P.Val)*100,2)
           fnFinalTab$Regulation <- ifelse(fnFinalTab$logFC <0,paste("UP_",fnConditions[1],"_DOWN_",fnConditions[2],collapse="",sep = ""),paste("DOWN_",fnConditions[1],"_UP_",fnConditions[2],collapse="",sep = ""))
           printOKMessage("      Results table .......................... OK")

           ####  Generando la grafica MD
           if(fnMDPlot)
           {
               regulationPoints<-data.frame(Exp=rep("Other",nrow(fnCountTable)),row.names=row.names(fnCountTable),stringsAsFactors = FALSE)
               deGenesUP = rownames(fnDeTab)[fnDeTab$adj.P.Val <= fnUmbral & fnDeTab$logFC <= -(fnUmbralFoldChange)]
               regulationPoints[deGenesUP,"Exp"]<-"UP"
               deGenesDOWN = rownames(fnDeTab)[fnDeTab$adj.P.Val <= fnUmbral & fnDeTab$logFC >= fnUmbralFoldChange]
               regulationPoints[deGenesDOWN,"Exp"]<-"DOWN"
               fnPlotFileName<-paste(fnFileName,"_plotMD",collapse="",sep = "")
               pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
               png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
               for(fnDevice in dev.list())
               {
                   dev.set(fnDevice)
                   plotMD(fnVoomedFitted, status=regulationPoints[,"Exp"],main=fnConditionsNames)
                   abline(h =c(-fnUmbralFoldChange, fnUmbralFoldChange), col = "blue") #c(-2, 2)
               }
               graphics.off()
               printOKMessage("      MD plot .......................... OK")
           }

           ####  Guardado de los datos en archivo
           if(TOP)
           {
               fnTopRows<-which((fnFinalTab$adj.P.Val<=fnUmbral) & (abs(fnFinalTab$logFC) >= fnUmbralFoldChange))
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
   else{
       printErrorMessage("      Limma (no replicates) .......................... Failed")
       return(0)
   }
}
