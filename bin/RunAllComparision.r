###  Autora:
###      Leticia Vega Alvarado.
###  Ultima actualizacion:
###      21/enero/2019
###  Descripcion:
###      Pipe line para realizar el calculo de Expresion diferencial por medio de tres metodos: edgeR, DESeq y NOISeq. El pipeline permite obtener el diagrama de Venn correspondiente
###      a la interseccion de los resultados de los 3 metodos.
###  Ejemplos:
### 21/09/2015 source("~/bin/DifferentialExp/RunAllComparision.r")
### 09/08/2016 RunAllComparisionMain("~/bin/DifferentialExp","~/Dropbox/Proyectos/Leonor_Perez/ExpDif_transcriptoma/Table_TF_All_Counts.txt","~/Dropbox/Proyectos/Leonor_Perez/ExpDif_transcriptoma",TOP=TRUE,fnUmbral=0.01,fnUmbralFoldChange=2, fnSelectMethods=12346)
### ("ColMS","ago42MS","ColMS","ago43MS","ColNaCl","ago42NaCl","ColNaCl","ago43NaCl")
### 10/08/2017 RunAllComparisionMain("~/bin/DifferentialExp","~/Temp/Flosan/modificado.txt","~/Temp/Flosan",TOP=TRUE,fnUmbral=0.01,fnUmbralFoldChange=2, fnSelectMethods=123456)
### RunAllComparisionMain("~/bin/DifferentialExp","~/Proyectos/Miguel/ConteosModificados.txt","~/Proyectos/Miguel",TOP=TRUE,fnUmbral=0.01,fnUmbralFoldChange=2, fnSelectMethods=13456)
### RunAllComparisionMain("~/bin/DifferentialExp","~/Proyectos/Babesia/Transc_All_Counts_reducido140218.txt","~/Proyectos/Babesia",TOP=TRUE,fnUmbral=0.05,fnUmbralFoldChange=1.5, fnSelectMethods=13456)

library(methods)
library(edgeR)

getPackageVersion<-function(fnMethods)
{
    fnMethods<-c(fnMethods,"ggplot2","UpSetR","corrplot","ComplexHeatmap")
    fnPackInfo<-installed.packages(fields=c("Package","Version"))
    fnPackInfo<-fnPackInfo[fnMethods,"Version"]
    fnPackageVersions<-paste(names(fnPackInfo),fnPackInfo,sep=" ")
    fnPackageVersions<-c(R.version.string,fnPackageVersions)
    return(fnPackageVersions)
}

loadSources<-function(fnMethods,fnProgamsPath)
{
    for(i in fnMethods){
        source(paste(fnProgamsPath,"/Run",i,".r",collapse="",sep = ""))
    }
    if("VennDiagram" %in% fnMethods){
        source(paste(fnProgamsPath,"/RunGetInfoFromFile.r",collapse="",sep = ""))
    }
    source(paste(fnProgamsPath,"/RunPrintMessage.r",collapse="",sep = ""))
    source(paste(fnProgamsPath,"/HeatmapPlots.r",collapse="",sep = ""))
}

createResultDir<-function(fnOutputPath,fnDEMethods)
{
    for(i in fnDEMethods){
        dir.create(paste(fnOutputPath,"/",i,"_Results",collapse="",sep = ""), showWarnings = FALSE, recursive = FALSE, mode = "0777")}
}

getCombinations<-function(fnSamplesName)
{
    fnCombVect<-NULL
    fnNumOfSamples=length(fnSamplesName)
    for (i in seq(1,fnNumOfSamples-1))
    {
        for(j in seq(i+1,fnNumOfSamples))
        {
            fnPair<-c(fnSamplesName[i],fnSamplesName[j])
            fnCombVect<-c(fnCombVect,fnPair)
        }
    }
    return(fnCombVect)
}

printRunInitInfo<-function(fnRunParameters,fnMethods)
{
    today <- Sys.Date()
    print("*************************  Running Date  *************************")
    print(today, format="%B %d %Y")
    print("*************************  Program call parameters  *************************")
    print(fnRunParameters)
    print("*************************  Software Versions  *************************")
    cat(paste(paste(getPackageVersion(fnMethods),collapse="\n"),"\n",sep="",collapse=""))
}

countFilterByCPM<-function(fnCountTable,fnFilterValue=1,fnConditions)
{
    fnSamplesName=factor(sub("_[a-zA-Z0-9]+$","",colnames(fnCountTable)))
    fnN<-min(table(fnSamplesName))
    ####  Filtrando los genes con pocas lecturas, en terminos del conteo por millon (cpm)
    fnFilterCountTable<-fnCountTable[rowSums(cpm(fnCountTable) >= fnFilterValue) >= fnN,]
    printOKMessage("      Filter rawcount table .......................... OK")
    return(fnFilterCountTable)
}

RunAllComparisionMain<-function(fnProgamsPath,fnInputFileName,fnOutputPath,fnCombinations=c(),TOP=FALSE,fnUmbral=0.01,fnFilterValue=1,fnUmbralFoldChange=1,fnSelectMethods=12345,fnFileGzipTar=TRUE,fnBatch=c())
{
   #####  Definicion de variables
   fnMethods<-c("edgeR","limma","NOISeq","DESeq2","DataAnalysis","VennDiagram")
   fnNomenclature<-as.integer(unlist(strsplit(as.character(fnSelectMethods),"")))
   fnDEMethods<-fnMethods[fnNomenclature[fnNomenclature>0 & fnNomenclature<5]]

   ##### Apertura del archivo .log
   sink(paste(fnOutputPath,"/","RunSummary.log",collapse="",sep = ""))
   
   ##### Imprimiendo en el archivo log la fecha de corrida, las versiones de paquetes utilizados y la línea de comandos con la que fue llamada la funcion RunAllComparisionMain
   fnMethodToPrint<-paste("RunAllComparisionMain(",fnProgamsPath,",",fnInputFileName,",",fnOutputPath,",fnCombinations=c(",paste(fnCombinations,collapse=",",sep=""),"),TOP=",TOP,",fnUmbral=",fnUmbral,",fnFilterValue=",fnFilterValue,",fnUmbralFoldChange=",fnUmbralFoldChange,",fnSelectMethods=",fnSelectMethods,",fnFileGzipTar=",fnFileGzipTar,",fnBatch=c(",paste(fnBatch,collapse=",",sep=""),"))",collapse="",sep="")
   printRunInitInfo(fnMethodToPrint,fnMethods[fnNomenclature[fnNomenclature!=5]])

   
   #####  Cargando los programas necesarios para el analisis y generando los directorios de salida
   loadSources(fnMethods[fnNomenclature],fnProgamsPath)
   createResultDir(fnOutputPath,fnMethods[fnNomenclature])

   #####  Lectura de la tabla de conteos y de las condiciones
   #if(!(is(fnCounts<-try(read.table(fnInputFileName,header=TRUE, row.names=1,sep="\t",quote=""),silent=T),"try-error")))
   if(class(fnCounts<-try(read.table(fnInputFileName,header=TRUE, row.names=1,sep="\t",quote=""),silent=T))!= "try-error")
   {
       print("*************************  Lectura de datos  *************************")
       printOKMessage("      Read count table .......................... OK")
       fnCounts<-fnCounts[ , order(sub("_[a-zA-Z0-9]+$","",colnames(fnCounts)))]
       fnAllCond=sub("_[a-zA-Z0-9]+$","",colnames(fnCounts))
       fnSamplesName=levels(factor(fnAllCond))
   
       if(5 %in% fnNomenclature){
           print("*************************  Running Data Analysis for all conditions  *************************")
           fnOutputFileName<-paste(fnOutputPath,"/DataAnalysis_Results/AllConditions",collapse="",sep = "")
           fnMethodToPrint<-paste("dataAnalysis(fnCounts,c(",paste(fnAllCond,collapse=",",sep=""),"),",fnOutputFileName,")",collapse="",sep="")
           print(fnMethodToPrint)
           dataAnalysis(fnCounts,fnAllCond,fnOutputFileName)
       }
   
       if(length(fnDEMethods)>0)
       {
           print("*************************  Running Differential Expresion  *************************")
          if(length(fnCombinations)==0){
              fnCombinations<-getCombinations(fnSamplesName)}
          fnPairCombination<-length(fnCombinations)
          for (i in seq(1,(fnPairCombination-1),by=2))
          {
              fnDataForVenn<-list()
              fnNamesVenn<-c()
              fnConditions<-which((fnAllCond == fnCombinations[i]) | (fnAllCond ==fnCombinations[i+1]))
              fnCounTable<-fnCounts[,fnConditions]
              fnCountTableFilter<-countFilterByCPM(fnCounTable,fnFilterValue=fnFilterValue,fnConditions=c(fnCombinations[i],fnCombinations[i+1]))
              for(k in 1:length(fnDEMethods))
              {
                  fnOutputPathMethod<-paste(fnOutputPath,"/",fnDEMethods[k],"_Results",collapse="",sep = "")
                  fnMethod<-paste("Run",fnDEMethods[k],"(fnCountTableFilter,fnOutputPathMethod,TOP=TOP,fnUmbral=fnUmbral,fnUmbralFoldChange=fnUmbralFoldChange, fnBatch=fnBatch,fnConditions=c(fnCombinations[i],fnCombinations[i+1]))",collapse="",sep = "")
                  fnVar<-try(eval(parse(text=fnMethod)),silent=TRUE)
                  if((class(fnVar)!="try-error"))
                  {
                      if(length(fnVar)>0)
                      {
                          fnDataForVenn[[fnDEMethods[k]]]<-fnVar
                          fnNamesVenn<-c(fnNamesVenn,fnDEMethods[k])
                      }
                  }
                  else{    printErrorMessage(paste(" ",fnDEMethods[k],"  .......................... Failed"),as.character(attr(fnVar,"condition")))}
              }
              if(6 %in% fnNomenclature){
                  fnVennDiagram<-try(VennDiag(fnDataForVenn,paste(fnOutputPath,"/VennDiagram_Results",collapse="",sep = ""),paste(fnCombinations[i],"vs",fnCombinations[i+1],collapse="",sep = "")),silent=TRUE)
                  if((class(fnVennDiagram)!="try-error"))
                  {
                      if(length(unlist(fnDataForVenn)>0))
                      {
                         if(is(try(RunGetInfoFromFile(fnOutputPath,fnNamesVenn,paste(fnCombinations[i],"vs",fnCombinations[i+1],collapse="",sep = "")),silent=TRUE),"try-error"))
                          {
                              printErrorMessage("      RunGetInfo .......................... Failed")
                          }
                          else{
                              printOKMessage("      RunGetInfo .......................... OK")
                          }
                          if(is(try(RunDataAbundance(fnOutputPath,fnNamesVenn,paste(fnCombinations[i],"vs",fnCombinations[i+1],collapse="",sep = ""),sort(unique(unlist(fnDataForVenn)))),silent=TRUE),"try-error"))
                          {
                              printErrorMessage("      RunDataAbundance .......................... Failed")
                          }
                          else{
                              printOKMessage("      RunDataAbundance .......................... OK")
                          }
                          if(is(try(RunDataLogFC(fnOutputPath,fnNamesVenn,paste(fnCombinations[i],"vs",fnCombinations[i+1],collapse="",sep = ""),sort(unique(unlist(fnDataForVenn)))),silent=TRUE),"try-error"))
                          {
                              printErrorMessage("      RunDataLogFC .......................... Failed")
                          }
                          else{
                              printOKMessage("      RunDataLogFC .......................... OK")
                          }
                          if(is(try(RunDataPval(fnOutputPath,fnNamesVenn,paste(fnCombinations[i],"vs",fnCombinations[i+1],collapse="",sep = ""),sort(unique(unlist(fnDataForVenn)))),silent=TRUE),"try-error"))
                          {
                              printErrorMessage("      RunDataPval .......................... Failed")
                          }
                          else{
                              printOKMessage("      RunDataPval .......................... OK")
                          }
                          
                          fnDEGenes<-paste(fnOutputPath,"/VennDiagram_Results/",fnCombinations[i],"vs",fnCombinations[i+1],"_Intesrsect_TOP_IDs.txt",collapse="",sep = "")
                          if(is(try(heatmapDEGenes(fnDEGenes,fnCountTableFilter,paste(fnOutputPath,"/VennDiagram_Results",collapse="",sep = ""),c(fnCombinations[i],fnCombinations[i+1]),"_Intersect"),silent=TRUE),"try-error"))
                          {
                              printErrorMessage("      RunHeatmap Intersect .......................... Failed")
                          }
                          else{
                              printOKMessage("      RunHeatmap Intersect .......................... OK")
                          }
                          fnDEGenes<-paste(fnOutputPath,"/VennDiagram_Results/",fnCombinations[i],"vs",fnCombinations[i+1],"_Union_TOP_IDs.txt",collapse="",sep = "")
                          if(is(try(heatmapDEGenes(fnDEGenes,fnCountTableFilter,paste(fnOutputPath,"/VennDiagram_Results",collapse="",sep = ""),c(fnCombinations[i],fnCombinations[i+1]),"_Union"),silent=TRUE),"try-error"))
                          {
                              printErrorMessage("      RunHeatmap Union .......................... Failed")
                          }
                          else{
                              printOKMessage("      RunHeatmap Union .......................... OK")
                          }
                      }
                  }
                  else{   printErrorMessage("      Venn Diagram .......................... Failed",as.character(attr(fnVennDiagram,"condition")))}
              }

          }
       }
       if(fnFileGzipTar)
       {
           setwd(fnOutputPath)
           fnFilestoCZip<-c(list.dirs(recursive=FALSE,full.names=FALSE),list.files(pattern="*.log"))
           tar("DiffExpAllResults.tar.gz",files=fnFilestoCZip,tar='tar',compression='gzip')
       }
   }
   else{   printErrorMessage("      Read count table .......................... Failed",as.character(attr(fnCounts,"condition")))}
   sink()
}

