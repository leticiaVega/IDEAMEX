RunNOISeq<- function(fnCountTable,fnOutputPath,TOP=FALSE,fnUmbral=0.01,fnUmbralFoldChange=1, fnBatch=c(), fnConditions)
{
   print("*************************  Running NOISeq  *************************")
   fnMethodToPrint<-paste("RunNOISeq(fnCounTable,",fnOutputPath,",TOP=",TOP,",fnUmbral=",fnUmbral,",fnUmbralFoldChange=",fnUmbralFoldChange,",fnConditions=c(",fnConditions[1],",",fnConditions[2],")",")",collapse="",sep="")
   print(fnMethodToPrint)
   
   ####  Cargando la biblioteca
   if(!("package:NOISeq" %in% search()))
   {
       library(NOISeq)
   }
   ####  Iicializacion de variables
   fnSamplesName=factor(sub("_[a-zA-Z0-9]+$","",colnames(fnCountTable)))
   fnExpDesign<-data.frame(fnSamples = colnames(fnCountTable),fnFactor =fnSamplesName)
   if(length(fnBatch))
   {
       fnExpDesign$Batch<-factor(fnBatch)
   }
   fnConditionsNames<-paste(fnConditions[1],"vs",fnConditions[2],collapse="",sep = "")
   fnFileName<-paste(fnOutputPath,"/",fnConditionsNames,collapse="",sep = "")
   fnPlotFileName<-paste(fnFileName,"_plotPCA",collapse="",sep = "")
   fnOutputFileName<-paste(fnFileName,".txt",collapse="",sep = "")
   fnOutputFileNameTop<-paste(fnFileName,"_TOP.txt",collapse="",sep = "")
   print("############")
   print(paste("Samples: ",fnConditionsNames))
   print("############")

   ####  Initializacion de un objecto de tipo NOISeq
   fnMyData<-try(NOISeq::readData( data=fnCountTable, factors=fnExpDesign),silent=T)
   if(!(is(fnMyData,"try-error")))
   {
       if(length(fnExpDesign$fnSamples)<=2)
       {
           fnMyResults <- noiseq(fnMyData, factor = "fnFactor",conditions=c(fnConditions[1],fnConditions[2]), k = NULL, norm = "n", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")
           printOKMessage("      Differential expression estimation.......................... OK")
       }
       else
       {
           if(length(fnBatch))
           {
               fnMyData=ARSyNseq(fnMyData,factor="Batch",batch=TRUE,norm="tmm",logtransf=FALSE)
               fnK=0
               fnLC=1
               fnNorm="n"
           }
           else{
               fnK=0.05
               fnLC=0
               fnNorm="tmm"
           }
           #### Graficando PCA
           fnMyPCA<-NOISeq::dat(fnMyData,type="PCA")
           pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
           png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
           for(fnDevice in dev.list())
           {
               dev.set(fnDevice)
               explo.plot(fnMyPCA, factor = "fnFactor")
           }
           graphics.off()
           printOKMessage("      PCA Plot .......................... OK")
           
           fnMyResults <-noiseqbio(fnMyData, k = 0.05, norm = "tmm", factor = "fnFactor",conditions=c(fnConditions[1],fnConditions[2]),lc = 0, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345,filter =0)
           printOKMessage("      Differential expression estimation.......................... OK")
       }
       ####  Obtencion de la tabla de resultados
       fnDETab<-fnMyResults@results[[1]]
   
       ####  Ordenando datos por nombre de genes para incorporar informacion adicional
       fnDETab <- fnDETab[order(row.names(fnDETab)),]

       ####  Calculo de la abundancia normalizada
       fnAbundance<-cbind(assayData(fnMyData)$exprs,tmm(assayData(fnMyData)$exprs, k=0.5, lc = 0))
       fnAbundance<-fnAbundance[order(row.names(fnAbundance)),]
       fnAbundanceFileName<-paste(fnFileName,"_Abundances.txt",collapse="",sep = "")
       write.table(tmm(assayData(fnMyData)$exprs), file=fnAbundanceFileName, sep="\t",quote=FALSE)
       printOKMessage("      Abundance estimation .......................... OK")
       
       fnFinalTab<-cbind(fnDETab[,],fnAbundance[,])
       fnFinalTab <- fnFinalTab[order(fnFinalTab$prob,decreasing=TRUE),]
       fnPvalFileName<-paste(fnFileName,"_pval.txt",collapse="",sep = "")
       write.table((1-fnFinalTab[,"prob",drop=FALSE]), file=fnPvalFileName, sep="\t",quote=FALSE)
   
       ####  Calculo de la certeza y la regulacion
       fnFinalTab$DE<-round(fnFinalTab$prob*100,2)
       if(length(fnExpDesign$fnSamples)>2)
       {
           fnFinalTab$Regulation <- ifelse(fnFinalTab$theta > 0, paste("UP_",fnConditions[1],"_DOWN_",fnConditions[2],collapse="",sep = ""),paste("DOWN_",fnConditions[1],"_UP_",fnConditions[2],collapse="",sep = ""))
           fnlogFCFileName<-paste(fnFileName,"_logFC.txt",collapse="",sep = "")
           write.table(fnFinalTab[,"log2FC",drop=FALSE], file=fnlogFCFileName, sep="\t",quote=FALSE)
       }
       else
       {
           fnFinalTab$Regulation <- ifelse(fnFinalTab$M > 0, paste("UP_",fnConditions[1],"_DOWN_",fnConditions[2],collapse="",sep = ""),paste("DOWN_",fnConditions[1],"_UP_",fnConditions[2],collapse="",sep = ""))
           fnlogFCFileName<-paste(fnFileName,"_logFC.txt",collapse="",sep = "")
           write.table(fnFinalTab[,"M",drop=FALSE], file=fnlogFCFileName, sep="\t",quote=FALSE)
       }
       printOKMessage("      Results table .......................... OK")
       
       ####  Generando la grafica DE
       fnPlotFileName<-paste(fnFileName,"_plotExpr",collapse="",sep = "")
       pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
       png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
       for(fnDevice in dev.list())
       {
           dev.set(fnDevice)
           DE.plot(fnMyResults, q=(1-fnUmbral), graphic = "expr", log.scale = TRUE)
       }
       graphics.off()
       printOKMessage("      DE plot .......................... OK")
       ####  Generando la grafica MD
       fnPlotFileName<-paste(fnFileName,"_plotMD",collapse="",sep = "")
       pdf(paste(fnPlotFileName,".pdf",collapse="",sep=""))
       png(paste(fnPlotFileName,".png",collapse="",sep=""),height=128, width=128,res=20)
       for(fnDevice in dev.list())
       {
          dev.set(fnDevice)
          DE.plot(fnMyResults, q=(1-fnUmbral), graphic = "MD")
       }
       graphics.off()
       printOKMessage("      MD plot .......................... OK")
       ####  Guardado de los datos en archivo
       if(TOP)
       {
           if(length(fnExpDesign$fnSamples)>2)
           {
               ##### Aqui hay que modificar tambien, porque no hay log2fc sino M
               fnTopRows<-which((fnFinalTab$prob>=(1-fnUmbral)) & (abs(fnFinalTab$log2FC) >= fnUmbralFoldChange))
           }
           else
           {
               fnTopRows<-which((fnFinalTab$prob>=(1-fnUmbral)) & (abs(fnFinalTab$M) >= fnUmbralFoldChange))
           }
           write.table(fnFinalTab[fnTopRows,],file=fnOutputFileNameTop, sep="\t",quote=FALSE)
           fnTopName<-row.names(fnFinalTab[fnTopRows,])
           printOKMessage("      Top file .......................... OK")
       }
       else
       {
           fnTopName=c()
       }
       write.table(fnFinalTab,file=fnOutputFileName, sep="\t",quote=FALSE)
       printOKMessage("      Final Results table to file .......................... OK")
       return(fnTopName)
   }
   else{
       printErrorMessage("      Objeto NOISeq .......................... Failed")
       return(fnMyData)
   }
}
