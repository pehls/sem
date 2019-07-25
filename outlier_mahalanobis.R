#######################################################################################
packages <- c('stringr','httpuv','lavaan','foreign','semTools','rJava','XLConnect','survey','tidyverse',
              'rdrop2','RCurl','matrixcalc','reshape2','graphics','mvoutlier','psych','heplots')
inst.pack <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
inst.pack(packages)


###### declara??o do modelo #####
model_completo <- '
  classroom=~INF_02_03+INF_02_04+INF_02_05+INF_02_06+INF_02_07
it_labs=~INF_03_01+INF_03_02+INF_03_03+INF_03_04+INF_03_06
specific_labs=~INF_04_01+INF_04_02+INF_04_03+INF_04_06
library=~INF_05_02+INF_05_04+INF_05_10+INF_05_09+INF_05_07+INF_05_21+INF_05_20
blackboard=~INF_06_01+INF_06_02+INF_06_04+INF_06_07+INF_06_08
program=~PROG_02_01+PROG_02_03+PROG_02_04+PROG_02_05
faculty=~PROF_02_01+PROF_02_02+PROF_02_03+PROF_02_04+PROF_02_05+PROF_02_06+PROF_02_07+PROF_02_08+PROF_02_15
disc_online=~ELEAR_02_01+ELEAR_02_02+ELEAR_02_03+ELEAR_02_04+ELEAR_02_05+ELEAR_02_06
coord=~COORD_02_01+COORD_02_02+COORD_02_03+COORD_02_04
financial_serv=~FSERV_01_01+FSERV_01_02+FSERV_01_03+FSERV_01_04
std_services=~SSERV_01_01+SSERV_01_02+SSERV_01_03+SSERV_01_04
image=~IMG_01_01+IMG_01_02+IMG_01_03
campus=~INF_01_01+INF_01_02+INF_01_04+INF_01_06+INF_01_09+INF_01_11+INF_01_13+INF_01_15+INF_01_12+INF_01_08
call=~CALL_01_01+CALL_01_02+CALL_01_03+CALL_01_04
employ=~EMPL_02_01+EMPL_02_02+EMPL_02_03+EMPL_02_04+EMPL_02_05+EMPL_02_06
inter=~INTL_01_01+INTL_01_02+INTL_01_03+INTL_01_04+INTL_01_05+INTL_01_07+INTL_01_06
infra=~campus+classroom+it_labs+specific_labs+library+blackboard
services=~financial_serv+std_services+call+employ+inter
#
SAT_00_02~infra+faculty+coord+program+disc_online+services+image
SAT_00_01~infra+faculty+coord+program+disc_online+services+image+SAT_00_02
NPS_01_01~SAT_00_01+image'
quebras <- c("Finished")
##########################
catch.lista.variaveis <- function(model_completo) {
  # agora vamos dividir entre constuctos e preditores
  pegar1 <- unlist(strsplit(model_completo,
                            split=c("\n")))
  pegar2 <- grep("=~", pegar1, value=T)
  pegar3 <- list()
  pegar4 <- c()
  z <- 1
  for(i in 1:length(pegar2)){
    
    pegar3[i] <- strsplit(pegar2[i],
                          split=c("=~"))
    pegar4[i] <- trimws(pegar2[[i]][1])
    pegar4    <- trimws(pegar4)
    
  }
  
  # separar os preditores dentro de cada constructo
  pegar5<-list()
  
  for(i in 1:length(pegar3)){
    
    ## vamos dividir em cada elemento da parte dos preditores
    pegar5[i] <- strsplit(pegar3[[i]][2],
                          split = c("\\+"))
    
  }
  
  # define a lista de variaveis do modelo, a partir dos constructos e das sozinhas
  lista_variaveis     <- unlist(pegar5)
  for (i in 1:length(lista_variaveis)) {
    lista_variaveis[i] <- trimws(lista_variaveis[i])
  }
  return(lista_variaveis) 
}

cleaning.db <- function(model_completo, var_lonely, dados) {
  lista_variaveis <- catch.lista.variaveis(model_completo)
  variaveis_do_modelo <- Reduce(intersect, list(lista_variaveis,colnames(dados)))
  variaveis_do_modelo <- c(variaveis_do_modelo,var_lonely)
  ## transformar as vaiaveis usadas em numericas
  indices_colunas <-  which(colnames(dados) %in% variaveis_do_modelo)
  remove_vars <- NULL
  for(i in 1:length(indices_colunas)){
    
    j <- indices_colunas[i]
    dados[,j] <-  as.numeric(as.character(dados[,j]))
    summary <- summary(dados[,j])
    if ((summary[7]/nrow(dados))>0.2) {
      remove_vars <- cbind(remove_vars, colnames(dados)[j])
    }
    
  }
  return(c(remove_vars))
  
}

do.remove<- function(model_completo_aux, variaveis_removidas) {
  ## primeiro passo eh dividir os models latentes
  model<-NULL
  model<-unlist(strsplit(model_completo_aux,
                         split=c("\n")))
  
  latent<-grep("=~", model, value=T)
  reg<-grep("~", model, value = T)
  reg<- reg[which(!reg %in%latent)]
  # agora vamos dividir entre constuctos e 'preditores'
  latentes<-list()
  nome_constructo<-c()
  
  for(i in 1:length(latent)){
    latentes[i]<-strsplit(latent[i],
                          split=c("=~"))
    nome_constructo[i]<-trimws(latentes[[i]][1])
    nome_constructo<-trimws(nome_constructo)
  }
  
  # print("dividir entre constuctos e preditores:ok")
  
  # separar os preditores dentro de cada constructo
  
  pred<-list()
  
  for(i in 1:length(latent)){
    
    ## vamos dividir em cada elemento da parte dos preditores 
    pred[i]<-strsplit(latentes[[i]][2],
                      split=c("\\+"))
    
  }
  varia<-(pred)
  ## colocar o nome do constructo em cada objeto
  names(varia)<-trimws(nome_constructo)
  
  ### VAMOS TIRAR AQUELAS VARIAVEIS QUE FORAM RETIRADAS
  ### NA INSPECAO DE MISSINGS
  k<-NULL
  var_aux<-NULL
  res<-NULL
  aux<-NULL
  
  for(i in 1:length(variaveis_removidas)){
    
    variaveis_PROBLEMA_MATRIZ_aux<-variaveis_removidas[i]
    res <- lapply(varia, function(ch) grep(variaveis_PROBLEMA_MATRIZ_aux, ch))
    aux<-sapply(res, function(x) length(x) > 0)
    
    if(sum(aux)>0){
      
      var_aux<-names(which(aux==T))
      k<-grep(variaveis_PROBLEMA_MATRIZ_aux, varia[[var_aux]])
      varia[[var_aux]][k]<-NA
      varia[[var_aux]]<-varia[[var_aux]][complete.cases(varia[[var_aux]])]
      
    }
    # variaveis_PROBLEMA_MATRIZ_aux<-NULL
    k<-NULL
    var_aux<-NULL
    res<-NULL
    aux<-NULL
  }
  
  HS.model<-NULL
  
  for(i in 1:length(varia)){
    
    va_lt<-NULL
    va_lt[i]<-paste(names(varia)[i],"=~")
    modelo<-NULL
    var_lat<-(paste(unlist(varia[[i]]),"+"))
    modelo<-paste(var_lat, collapse = " ")
    modelo<-substr(modelo, 1, nchar(modelo)-2)
    modelo<-paste(va_lt[i], modelo,"\n")
    HS.model<-paste(HS.model, modelo)
    
  }
  
  va_lt<-NULL
  constructos_excluidos <- NULL
  ## vamos tirar aquelas regressoes cujas variaveis independentes ficaram vazias
  HS.model_aux <- NULL
  HS.model_aux <-strsplit(HS.model,
                          split=c("\\n"))
  
  for(i in 1:length(HS.model_aux[[1]])){
    HS.model_aux2 <- NULL
    HS.model_aux2<-strsplit(trimws(HS.model_aux[[1]][i]),
                            split=c("\\~"))
    
    if(length(HS.model_aux2[[1]])<2){
      HS.model_aux[[1]][i]<-NA
      constructos_excluidos <- c(constructos_excluidos, sub(HS.model_aux2[[1]][1],pattern = "=",replacement = ""))
    }
  }
  HS.model_aux <- HS.model_aux[[1]][complete.cases(HS.model_aux)]
  for(i in 1:length(HS.model_aux)){
    HS.model_aux3 <- NULL
    HS.model_aux3<-strsplit(trimws(HS.model_aux[i]),
                            split=c("\\~"))
    var_lat<-((strsplit(HS.model_aux3[[1]][2],split = c("\\+"))))
    for (j in 1:length(var_lat[[1]])) {
      var_lat[[1]][j]<-trimws(var_lat[[1]][j])
    }
    for (k in 1:length(constructos_excluidos)) {
      constructos_excluidos[k] <- trimws(constructos_excluidos[k])
    }
    for (l in 1:length(constructos_excluidos)) {
      for (m in 1:length(var_lat[[1]])){
        if (!(is.na(var_lat[[1]][m])) && constructos_excluidos[l] == var_lat[[1]][m]) {
          var_lat[[1]][m] <- NA
        }
      }
    }
    HS.model_aux3[[1]][2] <- paste(var_lat[[1]][complete.cases(var_lat[[1]])], collapse=" + ")
    HS.model_aux[i]<-paste((HS.model_aux3[[1]]), collapse="~")
    
  }
  HS.model_aux <- HS.model_aux[complete.cases(HS.model_aux)]
  modelo2      <- NULL
  HS.model2    <- NULL
  
  for (i in 1:length(HS.model_aux)) {
    var_lat_aux2 <- NULL
    var_lat_aux2 <- (paste(unlist(HS.model_aux[i]),"\n"))
    modelo2      <- paste(var_lat_aux2, collapse = " ")
    HS.model2    <- paste(HS.model2, modelo2)
  }
  
  latent<-grep("=~", model, value=T)
  reg<-grep("~", model, value = T)
  reg<- reg[which(!reg %in%latent)]
  constructos.regressoes <- NULL
  for (i in 1:length(reg)) {
    constructos.regressoes <- str_split(reg,"~" )
  }
  var_lat <- NULL
  for (i in 1:length(constructos.regressoes)) {
    var_lat <- constructos.regressoes[[i]][2]
    var_lat <- str_split(var_lat, c("\\+"))
    for (j in 1:length(var_lat[[1]])) {
      var_lat[[1]][j] <- trimws(var_lat[[1]][j])
    }
    for (j in 1:length(constructos_excluidos)) {
      for (k in 1:length(var_lat[[1]])) {
        if (constructos_excluidos[j] == var_lat[[1]][k]) {
          var_lat[[1]][k] <- NA
        }
      }
      var_lat[[1]] <- var_lat[[1]][complete.cases(var_lat[[1]])]
    }
    var_lat <- paste(unlist(var_lat[[1]]), collapse="+")
    constructos.regressoes[[i]][2] <- var_lat
    reg[i] <- paste(unlist(constructos.regressoes[[i]]), collapse = "~")
  }
  HS.model2    <- paste(HS.model2, " # \n")
  
  HS.model2 <- paste(HS.model2,paste(unlist(reg), collapse = " \n "), sep="")
  
  return(HS.model2)
}
setwd("C:\\Users\\gabriel.pehls\\Dropbox\\UNIFACS2")
dbName <- "UNIFACS2_Satisfaction_2019.sav"
data_raw <<- read.spss(file          = dbName,
                       to.data.frame = T,
                       header        = T,
                       use.value.labels = F,
                       use.missings  = T)
#model_completo <- do.remove(model_completo, cleaning.db(model_completo, var_lonely = NULL, dados = data_raw))

model_completo <- '
  classroom=~INF_02_03+INF_02_04+INF_02_05+INF_02_06+INF_02_07
it_labs=~INF_03_01+INF_03_02+INF_03_03+INF_03_04+INF_03_06
specific_labs=~INF_04_01+INF_04_02+INF_04_03+INF_04_06
library=~INF_05_02+INF_05_04+INF_05_10+INF_05_09+INF_05_07+INF_05_21+INF_05_20
blackboard=~INF_06_01+INF_06_02+INF_06_04+INF_06_07+INF_06_08
program=~PROG_02_01+PROG_02_03+PROG_02_04+PROG_02_05
faculty=~PROF_02_01+PROF_02_02+PROF_02_03+PROF_02_04+PROF_02_05+PROF_02_06+PROF_02_07+PROF_02_08+PROF_02_15
coord=~COORD_02_01+COORD_02_02+COORD_02_03+COORD_02_04
financial_serv=~FSERV_01_01+FSERV_01_02+FSERV_01_03+FSERV_01_04
std_services=~SSERV_01_01+SSERV_01_02+SSERV_01_03+SSERV_01_04
image=~IMG_01_01+IMG_01_02+IMG_01_03
campus=~INF_01_01+INF_01_02+INF_01_04+INF_01_06+INF_01_09+INF_01_11+INF_01_13+INF_01_15+INF_01_12+INF_01_08
call=~CALL_01_01+CALL_01_02+CALL_01_03+CALL_01_04
employ=~EMPL_02_01+EMPL_02_02+EMPL_02_03+EMPL_02_04+EMPL_02_05+EMPL_02_06
inter=~INTL_01_01+INTL_01_02+INTL_01_03+INTL_01_04+INTL_01_05+INTL_01_07+INTL_01_06
infra=~campus+classroom+it_labs+specific_labs+library+blackboard
services=~financial_serv+std_services+call+employ+inter
#
SAT_00_02~infra+faculty+coord+program+services+image
SAT_00_01~infra+faculty+coord+program+services+image+SAT_00_02
NPS_01_01~SAT_00_01+image'
lista_variaveis     <- catch.lista.variaveis(model_completo)
variaveis_do_modelo <- Reduce(intersect, list(lista_variaveis,colnames(data_raw)))
data_raw<- as.data.frame(data_raw)
############## VFFF ####
detect.outliers <- function (data_raw, variaveis_do_modelo) {
  ##### imputar medias para missing values ####
  ## transformar as vaiaveis usadas em numericas
  indices_colunas <-  which(colnames(data_raw) %in% variaveis_do_modelo)
  
  for(i in 1:length(indices_colunas)){
    
    j <- indices_colunas[i]
    data_raw[,j] <-  as.numeric(as.character(data_raw[,j]))
    
  }
  
  indices_colunas <- variaveis_do_modelo
  
  # Finalmente, vamos adicionar a media nos NA's das variaveis
  # imputar medias gerais pra cada variavel
  data_raw$concatenar<-do.call(paste, as.list(data_raw[,quebras]))
  
  # percorrer cada individuo, em cada variavel, ver se esta nulo e colocar
  # a media da variavel conforme a vertical que pertence
  for(k in 1:length(indices_colunas)){
    
    j<-indices_colunas[k]
    
    teste_nan <- NULL
    teste_nan <- aggregate(data_raw[,j],
                           list(data_raw$concatenar), 
                           mean,
                           na.rm=T)
    
    svy.df <- NULL
    svy.df <- svydesign(id=~ResponseId, 
                        weights=~WEIGHT,
                        data=data_raw)
    
    media_geal_aux <- NULL
    media_geal_aux <- svymean(~data_raw[,j], 
                              design=svy.df,
                              na.rm=T,
                              deff=T)[1]
    
    teste_nan$x[is.nan(teste_nan$x)] <- media_geal_aux
    
    for(i in 1:nrow(data_raw)){
      
      vertical_aux<-data_raw$concatenar[i]
      
      if(is.na(data_raw[i,j])){
        #print(paste(i," ", j))
        data_raw[i,j] <- teste_nan$x[teste_nan$Group.1==vertical_aux]
        
      }
    }
  }
  
  data_raw$D <- sqrt(
    heplots::Mahalanobis(data_raw[,variaveis_do_modelo], 
                         colMeans(data_raw[,variaveis_do_modelo]), 
                         cov(data_raw[,variaveis_do_modelo]))
  )
  data_raw$outliers <- F
  data_raw$outliers[data_raw$D > sqrt(qchisq(0.999, df = length(variaveis_do_modelo)))] <- T
  data_raw$outlier_by_MD <- data_raw$outliers
  #table(data_raw$outliers)
  data_raw$var <- apply(data_raw[,variaveis_do_modelo], 1, var)
  data_raw$outliers[data_raw$var < 0.1] <- T
  data_raw$outlier_by_var<- F
  data_raw$outlier_by_var[data_raw$var < 0.1] <- T
  table(data_raw$outliers)
  #data_raw <- data_raw[which(data_raw$outliers == F),]
  x <- data_raw[,variaveis_do_modelo]
  require(mvoutlier)
  require(robustbase)
  covr <- covMcd(x, alpha=1)
  dist <- mahalanobis(x, center=covr$center, cov=covr$cov)
  s <- sort(dist, index=TRUE)
  if(ncol(x) > 2) {
    p <- princomp(x,covmat=covr)
    z <- p$scores[,1:2]
    sdprop <- (p$sd[1]+p$sd[2])/sum(p$sd)
    cat("Projection to the first and second robust principal components.\n")
    cat("Proportion of total variation (explained variance): ")
    cat(sdprop)
    cat("\n")
  }
  z <- as.data.frame(z)
  z$outliers <- data_raw$outliers
  z$outlier_byMD <- data_raw$outlier_by_MD
  z$outlier_byVAR <- data_raw$outlier_by_var
  write.csv2(z,"princcomp.csv")
  return(data_raw)
  
}

data_raw <- detect.outliers(data_raw, variaveis_do_modelo)

setwd("C:\\Users\\gabriel.pehls\\Google Drive\\TCC\\Apresentação")
write.csv2(data_raw, "dados_outliers.csv")
