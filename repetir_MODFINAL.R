rm(list=ls())

########################################################################################
###| Complex sampling analysis of SEM models--------|###################################
###| Daniel Oberski, 2015-11-03---------------------|###################################
###| Updated by Pehls, 2017-12-08-------------------|###################################
########################################################################################
lavaan.survey <- 
  function(lavaan.fit, survey.design, 
           estimator=c("MLM", "MLMV", "MLMVS", "WLS", "DWLS", "ML"),
           estimator.gamma=c("default","Yuan-Bentler")) {
    
    # Not all estimators in lavaan make sense to use here, therefore matching args
    #print("inicializando estimator")
    estimator <- match.arg(estimator) 
    if(estimator=="ML") warning("Estimator 'ML' will not correct standard errors and chi-square statistic.")
    #print("inicializando teste estimator")
    estimator.gamma <- match.arg(estimator.gamma) # Smoothing Gamma or not
    #print("inicializando teste lavaan names")
    # Names of the observed variables (same for each group)
    ov.names <- lavaanNames(lavaan.fit, type="ov", group=1)
    #print("inicializando dplus (matrix duplication)")
    # The MP-inverse duplication matrix is handy for removing redundancy
    Dplus <- lavaan::lav_matrix_duplication_ginv( length(ov.names) )
    #print("tamanho do Dplus:")
    #print(object.size(Dplus))
    #print("inicializando criacao formula")
    # Create a formula that includes all observed variables for svymean
    ov.formula <- as.formula(paste("~", paste(ov.names, collapse="+")))
    #print("inicializando inspecao ngroups")
    # <no. group>-sized lists that will contain the asy. covariance matrix,
    #  and sample covariance matrix and mean vector
    ngroups <- lavInspect(lavaan.fit, "ngroups")
    #print("inicializando Gamma vector")
    Gamma <- vector("list", ngroups)
    #p#rint("tamanho do Gamma vector:")
    #print(object.size(Gamma))
    #print("inicializando sample.cov")
    sample.cov <- vector("list", ngroups)
    #print("inicializando sample.mean")
    sample.mean <- vector("list", ngroups)
    #print("inicializando sample.nobs")
    sample.nobs <- vector("list", ngroups)
    
    #print("inicializando for ngroups")
    for(g in seq(ngroups)) {
      if(ngroups > 1) {
        # Use survey::subset to create data groups
        #print("inicializando survey.design.g")
        survey.design.g <- 
          subset(survey.design, eval(parse(text=sprintf("%s == '%s'", 
                                                        lavInspect(lavaan.fit, "group"), 
                                                        lavInspect(lavaan.fit, "group.label")[[g]]))))
      } 
      else { # In case of no groups, just use the original survey design object.
        #print("inicializando survey.design.g no groups")
        survey.design.g <- survey.design  
      }
      
      # Function that takes survey design and returns the Gamma & observed moments
      get.stats.design <- function(survey.design.g, sample.nobs.g) {
        #print("inicializando matrix transformation")
        ##essa operacao e demorada
        sample.cov.g <- as.matrix(svyvar(ov.formula, design=survey.design.g, na.rm=TRUE))  
        #print("inicializando covariances matrix as attribute")
        # survey package returns the variance matrix of the (co)variances as attr:
        Gamma.cov.g <- attr(sample.cov.g, "var")
        #print("inicializando remocao covariances")
        ##operacao pode ser demorada
        # Remove (co)variances wrt redundant elements of covs; not used by lavaan. 
        Gamma.cov.g <- Dplus %*% Gamma.cov.g %*% t(Dplus)
        #print("inicializando remocao mean vector")
        # Same for mean vector
        sample.mean.g <- svymean(ov.formula, design=survey.design.g, na.rm=TRUE)  
        Gamma.mean.g <- attr(sample.mean.g, "var")
        # Join asy. variance matrices for means and covariances
        # TODO add offdiag
        #print("inicializando lavaan matrix")
        Gamma.g <- lavaan::lav_matrix_bdiag(Gamma.mean.g, Gamma.cov.g)
        Gamma.g <- Gamma.g * sample.nobs.g # lavaan wants nobs * Gamma.
        # Since the above nonparametric estimate of Gamma can be unstable, Yuan
        # and Bentler suggested a model-smoothed estimate of it, optional here:
        if(estimator.gamma == "Yuan-Bentler") {
          #print("inicializando yuan-bentler residuals")
          r <- get.residuals(lavaan.fit) # Iff these asy = 0, all will be well...
          Gamma.g <- Gamma.g + (sample.nobs.g/(sample.nobs.g - 1)) * (r %*% t(r))
        }
        # Get rid of attribute, preventing errors caused by lazy evaluation
        # (This has to be at the end or lazy evaluation mayhem will ensue)
        attr(sample.cov.g, "var") <- NULL
        tmp  <- as.vector(sample.mean.g)
        names(tmp) <- names(sample.mean.g)
        sample.mean.g <- tmp
        #print("inicializando list gamma")
        list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
      }
      # The data may be a list of multiply imputed datasets
      if(!any(class(survey.design.g) == "svyimputationList")) {
        #print("inicializando lavInspect lavaan.fit")
        # If no imputations, just use usual no. observations and asy variance
        sample.nobs.g <- lavInspect(lavaan.fit, "nobs")[[g]] 
        #print("tamanho do sample.nobs.g:")
        #print(object.size(sample.nobs.g))
        #print("inicializando get.stats")
        stats <- get.stats.design(survey.design.g, sample.nobs.g)
        #print("tamanho do stats:")
        #print(object.size(stats))
        
      } 
      else { # In case of multiply imputed data
        # Not only can nobs differ from lavaan.fit, but also per imputation
        sample.nobs.g <- get.sample.nobs(survey.design.g, lavInspect(lavaan.fit, "group"))
        
        # Retrieve point and variance estimates per imputation, 
        #    [TODO: this line will not be correct when nobs differs over imputations]
        stats.list <- lapply(survey.design.g[[1]], get.stats.design, sample.nobs=sample.nobs.g)
        m  <- length(stats.list) # no. imputation
        
        # Point estimates are average over imputations
        sample.cov.list <- lapply(stats.list, `[[`, 'sample.cov.g')
        sample.cov.g <- Reduce(`+`, sample.cov.list) / m
        cov.df <- Reduce(`rbind`, lapply(sample.cov.list, lavaan::lav_matrix_vech))
        sample.mean.list <- lapply(stats.list, `[[`, 'sample.mean.g')
        sample.mean.g <- Reduce(`+`, sample.mean.list) / m
        mean.df <- Reduce(`rbind`, sample.mean.list)
        #print("inicializando gamma reduction/calcs")
        # Variance estimates depend on within- and between-imputation variance:
        Gamma.within  <- Reduce(`+`, lapply(stats.list, `[[`, 'Gamma.g')) / m
        Gamma.between <- cov(cbind(mean.df, cov.df))
        Gamma.g <- Gamma.within + ((m + 1)/m) * Gamma.between
        
        # set stats with multiple imputation point and variance estimates
        stats <- list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
      }
      #print("inicializando juncao lists")
      # Augment the list for this group
      Gamma[[g]] <- stats$Gamma.g
      sample.cov[[g]] <- stats$sample.cov.g
      sample.mean[[g]] <- stats$sample.mean.g
      sample.nobs[[g]] <- sample.nobs.g
      stats <- NULL
    } # End of loop over groups
    #print("end of loop ngroup")
    Dplus <- NULL
    #print("preenchendo new.call no lavaan.survey")
    new.call <- lavInspect(lavaan.fit, "call")
    new.call$data <- NULL                # Remove any data argument
    new.call$sample.cov <- sample.cov    # Set survey covariances
    new.call$sample.mean <- sample.mean  # Set survey means
    new.call$sample.nobs <- sample.nobs  
    new.call$estimator <- estimator  # Always use Satorra-Bentler or WLS estimator
    #print("tamanho do new.call")
    #print(object.size(new.call))
    if(substr(estimator, 1, 2) == "ML") { # ML, robust options
      # Set asymptotic covariance matrix of sample means and covariances
      new.call$NACOV <- Gamma  
    }
    if(estimator %in% c("WLS", "DWLS")) {
      # Weighted Least Squares, adjust the weight matrix: MP inverse of Gamma
      # Note that Gamma may be singular.
      new.call$WLS.V <- lapply(Gamma, ginv)
    }
    #print("preenchendo new.fit")
    ##outra operacao demorada
    if (is.positive.definite((cov(new.call$sample.cov[[1]])))) { # in a case, the new.call wasn't positive definite, and lavaan can't run, so
      # if the new.call was positive definite
      # print("tamanho do lavaan.fit:")
      # print(object.size(lavaan.fit))
      # print("tamanho do estimador Gamma:")
      #  print(object.size(Gamma))
      new.fit <- eval(as.call(new.call)) # Run lavaan with the new arguments
    } else {
      return(lavaan.fit) # the new.fit was the previous fit
    }
    new.call <- NULL
    if(estimator %in% c("WLS", "DWLS")) return(new.fit) # We are done for WLS
    
    # For ML with robust se's, check that a possibly singular Gamma has not
    # created dependencies in the parameter estimates.
    # (Code below should really be implemented in lavaan...)
    evs.too.small <- sapply(Gamma, function(Gamma.g) {
      any(eigen(Gamma.g, only.values=TRUE)$values < (.Machine$double.eps*100))
    })
    Gamma <- NULL
    if(any(evs.too.small)) {
      V.est <- lavaan::vcov(new.fit)
      if(any(Re(eigen(V.est, only.values=TRUE)$values) < (.Machine$double.eps*100))) {
        long.string  <- sprintf("Some of the standard errors may not be trustworthy.
                                Some of the observed covariances or means are
                                collinear, and this has generated collinearity in your
                                parameter estimates.  This may be a sample size issue,
                                missing data problem, or due to having too few
                                clusters relative to the number of parameters. Problem
                                encountered in group(s) %s",
                                paste(which(evs.too.small), collapse=", "))
        
        warning(strwrap(long.string, width=9999, simplify=TRUE))#gotta love it
      }
    }
    #print("finalizando new.fit")
    # print("tamanho do new.fit final:")
    #  print(object.size(new.fit))
    new.fit
  }

# Obtain residuals from a lavaan fit object, concatenating means w/ covariances
# (used in Yuan-Bentler correction)
get.residuals <- function(fit) {
  r  <- lavaan::residuals(fit)
  c(r$mean, lavaan::lav_matrix_vech(r$cov))
}

# Obtain sample size from multiply imputed svydesign object.
# In case sample size differs over imputations, takes median over imputations.
# TODO: Does not work with multiple group yet.
get.sample.nobs  <- function(svy.imp.design, group=NULL) {
  nobs.imp <- lapply(svy.imp.design[[1]], function(des) {nrow(des$variables)})
  return(median(unlist(nobs.imp)))
}

# Use the pFsum function from the survey package to obtain p value for the 
#   overall model fit using an F reference distribution where the 
#   denominator degrees of freedom is the design degrees of freedom.  
# An anonymous reviewer for J Stat Software suggested that 
#  "in surveys with small numbers of primary sampling units this sort of 
#   correction has often improved the 
#   behavior of tests in other contexts."
# The eigenvalues of the U.Gamma matrix will be the coefficients in the 
#   mixture of F's distribution (Skinner, Holt & Smith, pp. 86-87).
pval.pFsum <- function(lavaan.fit, survey.design, method = "saddlepoint") {
  # Check that Satorra-Bentler or Satterthwaite adjustment is present
  if(!lavInspect(lavaan.fit, "options")$test %in% 
     c("satorra.bentler", "mean.var.adjusted", "Satterthwaite")) {
    stop("Please refit the model with Satorra-Bentler (MLM) or Satterthwaite (MLMVS) adjustment.") 
  }
  
  UGamma <- lavTech(lavaan.fit, "ugamma")
  real.eigen.values <- Re(eigen(lavTech(lavaan.fit, "ugamma"), only.values = TRUE)$values)
  
  return(survey::pFsum(x=fitMeasures(lavaan.fit, "chisq"), df=rep(1, length(real.eigen.values)), 
                       a=real.eigen.values, ddf=survey::degf(survey.design), lower.tail=FALSE,
                       method=method))
}
################################ End of function #######################################

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
    if (!is.null(constructos_excluidos)) {
      for (l in 1:length(constructos_excluidos)) {
        for (m in 1:length(var_lat[[1]])){
          if (!(is.na(var_lat[[1]][m])) && constructos_excluidos[l] == var_lat[[1]][m]) {
            var_lat[[1]][m] <- NA
          }
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
  if (length(constructos.regressoes) > 0) {
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
  }
 if (length(reg)>0) {
   HS.model2    <- paste(HS.model2, " # \n")
   
   HS.model2 <- paste(HS.model2,paste(unlist(reg), collapse = " \n "), sep="")
 }
 
  
  return(HS.model2)
}




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
  #table(data_raw$outliers)
  data_raw$var <- apply(data_raw[,variaveis_do_modelo], 1, var)
  data_raw$outliers[data_raw$var < 0.1] <- T
  #table(data_raw$outliers)
  data_raw <- data_raw[which(data_raw$outliers == F),]
  return(data_raw)
  
}



do.boot <- function(IES,
                    quebras, n_quebra, pular,
                    model_completo, constructs, var_dep_eq, segunda_ordem, var_lonely, cod_dropbox) {
  packages <- c('boot','foreign','doParallel','stringr','httpuv',
                'lavaan','semTools','survey',
                'tidyverse','rdrop2','RCurl','matrixcalc',
                'reshape2','graphics','mvoutlier','psych','heplots')
  inst.pack <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  if (Sys.getenv("JAVA_HOME")!="")
    Sys.setenv(JAVA_HOME="")
  library(rJava)
  library(XLConnect)
  inst.pack(packages)
  
  ## catch year for name's
  year <- 2019
  #TROCAR AQUI CASO OFFLINE
  
  Estudo <- paste(IES, paste("Satisfaction", year, sep="_"), sep="_")
  dbName <- paste(Estudo, ".sav", sep="")
  
  ## Download from dropbox
  dl_from_dropbox <- function(x, key) {
    require(RCurl)
    bin <- getBinaryURL(paste0("https://dl.dropboxusercontent.com/s/", key, "/", x),
                        ssl.verifypeer = FALSE)
    con <- file(x, open = "wb")
    writeBin(bin, con)
    close(con)
    message(noquote(paste(x, "read into", getwd())))
  }
   # dl_from_dropbox(dbName,cod_dropbox)
  
  data_raw <<- read.spss(file          = dbName,
                         to.data.frame = T,
                         header        = T,
                         use.value.labels = F,
                         use.missings  = T)
  #pegar as variable labels
  var_labels <- data.frame(attr(data_raw, 'variable.labels'))    %>%
    cbind(., rownames(.) )
  
  rownames(var_labels) <- NULL
  
  colnames(var_labels) <- c("label","codigo")
  data_raw_aux_t2b <<- data_raw
  # baixar os dados com labels e filtrar conforme a quebra
  data_value_labels <<- read.spss(file          = dbName,
                                  to.data.frame = T,
                                  header        = T,
                                  use.value.labels = T,
                                  use.missings  = T)    
  Estudo_aux <- Estudo
  for(q in n_quebra:length(quebras)){
    #for(q in 1:1){
    ## quebra 'q'
    #q=1
    quebra  <- quebras[q]
    niv_aux <- droplevels(data_value_labels[,quebras[q]])
    niv     <- levels(niv_aux)
    Estudo <- paste(Estudo_aux, quebra, sep="_")
    # monitora tempo
    # tempo_inico_proc_quebra[q] <- proc.time()
    n <- 1
    if (q==n_quebra) {
      n <- pular
    }
    
    for(lev in n:length(niv)){
      #lev=1
      ## Identifica o nivel 'niv' da quebra 'q'
      nivel_aux <- droplevels(as.factor(data_raw[,quebras[q]]))
      nivel     <- levels(nivel_aux)[lev]
      
      # define a lista de variaveis do modelo, a partir dos constructos e das sozinhas
      lista_variaveis     <- catch.lista.variaveis(model_completo)
      variaveis_do_modelo <- Reduce(intersect, list(lista_variaveis,colnames(data_raw)))
      variaveis_do_modelo <<- c(variaveis_do_modelo,var_lonely)
      
      ## transformar as vaiaveis usadas em numericas
      indices_colunas <-  which(colnames(data_raw) %in% variaveis_do_modelo)
      
      ## Identifica o nome do nivel da quebra
      nivel_label_aux <- droplevels(as.factor(data_value_labels[,quebras[q]]))
      nivel_label     <- levels(nivel_label_aux)[lev]
      print("Iniciando o nivel ")
      print(nivel_label)
      ## verifica se o nivel ? valido para a analise e se a base nao ? muito pequena (n < 30) ou o numero de obs ? menor que o de variaveis do modelo
      if (is.na(nivel) ) break
      if ( (count(data_raw[data_raw[,quebra]==nivel,]) < 31 ) || ((count(data_raw[data_raw[,quebra]==nivel,]) < (length(indices_colunas)*10) ))){
        
        tamanho_base <- count(data_raw[data_raw[,quebra]==nivel,])
        print(" TAMANHO DA BASE INSUFICIENTE PARA ANaLISE ")
        print(paste(" n = ",tamanho_base[[1]]))
        
        fileName <- paste(Estudo,
                          paste(nivel_label,
                                paste(Sys.Date(),"witherror.xlsx"),sep="_"),sep="_")  %>%
          stringr::str_replace_all(.,"/","")  %>%
          stringr::str_replace_all(.," ","")  %>%
          stringr::str_replace_all(.,"-","")
        fileXls <- fileName
        unlink(fileXls, recursive = FALSE, force = FALSE)
        exc <- loadWorkbook(fileXls, create = TRUE)
        
        createSheet(exc,'Observacoes')
        saveWorkbook(exc)
        input <- "QUEBRA TORNA TAMANHO DA BASE INSUFICIENTE PARA ANaLISE"
        writeWorksheet(exc, tamanho_base, sheet ='Observacoes',header=T,rownames = T, startRow = 1, startCol = 2)
        writeWorksheet(exc, input, sheet ='Observacoes',header=T,rownames = T, startRow = 1, startCol = 3)
        saveWorkbook(exc)
        drop_upload(fileName, path = IES)
        
      }else {
        z <- 1
        data_aux_top <- NULL
        ##### CRIAcaO DE BANCOE MODELO AUXILIAR PARA ANaLISE DE QUEBRAS #####
        # cira banco auxiliar para verificacao de vazias no top2box
        variaveis_do_modelo <<- Reduce(intersect, list(lista_variaveis,colnames(data_raw)))  %>% 
          c(.,var_lonely)
        data_aux_top        <<- data_raw_aux_t2b[data_raw_aux_t2b[,quebra]==nivel,]          %>%
          select(variaveis_do_modelo, "WEIGHT", "ResponseId", quebras)
        
        model_completo_aux <- model_completo
        
        variaveis_removidas <- NULL
        
        variaveis_removidas <- c(variaveis_removidas)
        
        if (length(variaveis_removidas) == 0){
          mapeamente_vars_removidas <- cbind(variaveis= "Nenhuma",motivos="Variavel_Vazia", momento=z)
        }else{
          mapeamente_vars_removidas <- cbind(variaveis= variaveis_removidas,motivos="Variavel_Vazia/mais de 20% de missings", momento=z)
          model_completo_aux   <- model_completo
          
          model_completo_aux <- do.remove(model_completo_aux, variaveis_removidas)
        }
        
        # redefine varaveis do modelo excluindo variaveis nulas
        
        
        lista_variaveis <- catch.lista.variaveis(model_completo_aux)
        
        variaveis_do_modelo <- Reduce(intersect, list(lista_variaveis,colnames(data_raw)))
        variaveis_do_modelo <- c(variaveis_do_modelo,var_lonely)
        
        
        # cria banco auxiliar com a quebra para o raw data e para o top2box
        dados  <- NULL
        dados  <- data_raw[which(data_raw[,quebra]==nivel),]
        dados  <- dados[,c(var_dep_eq, "ResponseId",
                           variaveis_do_modelo,"WEIGHT", quebras)]
        dados <- as.data.frame(dados)
        
        
        dados <-detect.outliers(data_raw, variaveis_do_modelo)
        data_raw <- dados
        
        #registerDoParallel(cores = detectCores())
       
        
        for(i in 1:5) #%dopar%
          do.sem(indice = i, dados = dados,
                 var_labels = var_labels,
                 nivel = nivel,
                 nivel_label = nivel_label,
                 model_completo_aux = model_completo_aux,
                 mapeamente_vars_removidas = mapeamente_vars_removidas,
                 constructs = constructs,
                 var_dep_eq = var_dep_eq,
                 segunda_ordem = segunda_ordem,
                 var_lonely = var_lonely,
                 IES = IES)
      # 
     
   

  
  }
      var_retirada_porT2B <- NULL
      
      resumo_modelo = NULL
      
      r2_bancofim_v1 = NULL
      
      modelo_final_final = NULL
      
      est_modelo_impacto = NULL
      
      var_fora = NULL
      
      modelo_partial = NULL
      
      est_modelo_final_reg_val = NULL
      
      input = NULL
    }
  }
}


do.sem <- function(indice, dados, var_labels,
                   nivel, nivel_label, 
                   model_completo_aux, mapeamente_vars_removidas,
                   constructs, var_dep_eq, segunda_ordem, var_lonely, 
                   IES) {
  print(paste("iniciando paralelamente vers?o", indice, sep="  "))
  ## catch year for name's
  data <- date()
  year <- 2019
  hour <- (strsplit(data," ")[[1]])[4]
  
  Estudo <- paste(IES, paste("Satisfaction", year, sep="_"), sep="_")
  dbName <- paste(Estudo, ".sav", sep="")
  Estudo <- paste(Estudo,hour, sep='_') %>% stringr::str_replace_all(.,":","_")
  
  ## load/install all packages
  
  packages <- c('stringr','httpuv','lavaan','foreign','semTools','survey','tidyverse',
                'rdrop2','RCurl','matrixcalc','reshape2')
  inst.pack <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  if (Sys.getenv("JAVA_HOME")!="")
    Sys.setenv(JAVA_HOME="")
  library(rJava)
  library(XLConnect)
  inst.pack(packages)
  
  ## Download from dropbox
  dl_from_dropbox <- function(x, key) {
    require(RCurl)
    bin <- getBinaryURL(paste0("https://dl.dropboxusercontent.com/s/", key, "/", x),
                        ssl.verifypeer = FALSE)
    con <- file(x, open = "wb")
    writeBin(bin, con)
    close(con)
    message(noquote(paste(x, "read into", getwd())))
  }
  
  # dl_from_dropbox(".httr-oauth","1z9mwuxpqjsw8vl")
  # dl_from_dropbox("token.rds","nl74tpnq9krj54z")
  token <- readRDS("token.rds")
  drop_acc(dtoken = token)
  # como determinar o tamanho da base
  lista_variaveis <- catch.lista.variaveis(model_completo_aux)
  
  variaveis_do_modelo <- Reduce(intersect, list(lista_variaveis,colnames(dados)))
  variaveis_do_modelo <- c(variaveis_do_modelo,var_lonely)
  dados <-detect.outliers(data_raw, variaveis_do_modelo)
  data_raw <- dados
  tamanho <- 500
  data_raw <- (dados[sample(x = 1:nrow(dados), size = tamanho),])
  dados <- data_raw
  # baixar banco de dados sem labels e filtrar conforme a quebra ####
  
  data_raw_aux_t2b <<- data_raw
  
  
  svy.df  <-  svydesign(id=~ResponseId, 
                        weights=~WEIGHT,
                        data=data_raw)
  
  
  
  ## Vamos imputar valores para os missing values
  nome_constructo <- 0
  valat_originais <- 0
  diferentes      <- 0
  
  # define a lista de variaveis do modelo, a partir dos constructos e das sozinhas
  lista_variaveis     <- catch.lista.variaveis(model_completo_aux)
  variaveis_do_modelo <- Reduce(intersect, list(lista_variaveis,colnames(data_raw)))
  variaveis_do_modelo <<- c(variaveis_do_modelo,var_lonely)
  
  ## transformar as vaiaveis usadas em numericas
  indices_colunas <-  which(colnames(data_raw) %in% variaveis_do_modelo)
  
  for(i in 1:length(indices_colunas)){
    
    j <- indices_colunas[i]
    data_raw[,j] <-  as.numeric(as.character(data_raw[,j]))
    
  }
  
  ########################################################################################
  ############              ANTES: CALCULAR O T2B PARA CADA VARIAVEL          ############
  ########################################################################################
  
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
  
  # Atualizar a ponderacao os resultados
  svy.df<-svydesign(id=~ResponseId, 
                    weights=~WEIGHT,
                    data=data_raw)
  
  variancia_geral_nome <- NULL
  variancia_geral_aux   <- NULL
  
  for(k in 1:length(indices_colunas)){
    
    j<-indices_colunas[k]
    
    variancia_geral_nome[k] <- indices_colunas[k]
    variancia_geral_aux[k]   <- svyvar(~data_raw[,j], 
                                       design=svy.df,
                                       na.rm=T,
                                       deff=T)[1]
    #print(variancia_geral_aux[k])
  }
  
  variancia_geral <- cbind(indices_colunas, variancia_geral_aux)
  
  ############                     FIM DA SEcaO                               ############
  ########################################################################################
  
  ########################################################################################
  ############           CaLCULO DA ANaLISE FATORIAL E REGRESS?ES             ############
  ############                              AUTOMaTICO                        ############
  ############ nao sao necessarias modificac?es                               ############
  ########################################################################################
  
  # A gente vai fazer para cada quebra necessaria
  
  tempo_inicial_total <- proc.time()
  tempo_inico_proc_quebra    <- NULL
  tempo_final_proc_quebra    <- NULL
  Estudo_aux <- Estudo
  dados <- data_raw
  n_obs <- nrow(dados)
  
  data_aux_top <<- dados %>%
    select(variaveis_do_modelo, "WEIGHT", "ResponseId", quebras)
  # executado na linha 136
  for(i in 1:length(variaveis_do_modelo)){
    
    data_aux_top[,i]<- ifelse(data_aux_top[,i]>3,"T2B", 
                              ifelse(data_aux_top[,i]<3,"B2B","Neutral"))
    
  }
  
  print(paste("Equac?es Estruturais para Quebra",nivel_label))
  
  svy.df <- NULL
  svy.df <- svydesign(id=~ResponseId, 
                      weights=~WEIGHT,
                      data=dados)
  data_aux_top <<- data_raw_aux_t2b %>%
    select(variaveis_do_modelo, "WEIGHT", "ResponseId", quebras)
  # executado na linha 136
  for(i in 1:length(variaveis_do_modelo)){
    
    data_aux_top[,i]<- ifelse(data_aux_top[,i]>3,"T2B", 
                              ifelse(data_aux_top[,i]<3,"B2B","Neutral"))
    
  }
  
  ############################################################################
  ######### BANCO AUXILIAR PARA O TOP 2 BOX      #############################
  ############################################################################
  
  i<-1
  # Ponderar os resultados
  svy.df.auxtop<-NULL
  svy.df.auxtop<-svydesign(id=~ResponseId, 
                           weights=~WEIGHT,
                           data=data_aux_top)
  df.top2box<-data.frame()
  i<<-1
  for(i in 1:length(variaveis_do_modelo)){
    nome_aux_t2b <- NULL    
    nome_aux_t2b <- colnames(data_aux_top[i])
    aux_t2box    <- NULL
    aux_t2box    <- as.data.frame(prop.table(svytable(~eval(as.symbol(colnames(data_aux_top)[i])),
                                                      design = svy.df.auxtop))*100)
    aux_ver_exist<-NULL
    aux_ver_exist<-('T2B' %in% aux_t2box[,1])
    
    if(aux_ver_exist==FALSE){
      aux_t2box[,1]                 <- as.numeric(aux_t2box[,1])
      aux_t2box[nrow(aux_t2box)+1,1]<- 'T2B'
      aux_t2box[nrow(aux_t2box)+1,2]<- 0
    }
    
    df.top2box_aux<-NULL
    df.top2box_aux<-cbind(nome_aux_t2b,aux_t2box$Freq[aux_t2box[,1]=="T2B"])
    df.top2box    <-rbind(df.top2box,df.top2box_aux)
  }
  
  colnames(df.top2box) <-c("variaveis", "T2B_result")
  df.top2box$T2B_result<-as.numeric(as.character(df.top2box$T2B_result))
  df.top2box[which(is.na(df.top2box$T2B_result)),2]<-0
  
  ###############################################################################
  #### o simbolo "=~" significa que a a variavel latente y se manifesta pelas variaveis
  #### observadas x1, x2,...xk
  soh_val <- NULL
  soh_val <- sapply(strsplit(model_completo_aux, split= "\\#"),
                    function(x) x[length(x)])
  
  HS.model <- NULL
  HS.model <- substr(model_completo_aux,1,
                     regexpr("#",model_completo_aux)-1)
  
  
  #########################################################################
  ### VAMOS TESTAR SE A MATRIZ DE COVARIANCIA EH DEFINIDA POSITIVA #######
  # print("Executando teste de matriz de covariancia definida positiva...")
  teste_covariancia <- (is.positive.definite(cov(dados[,variaveis_do_modelo])) 
                        & n_obs>length(variaveis_do_modelo))
  
  #### Se for FALSE,  entao a matriz de covariancia nao eh definida, positiva.
  #### E a modelagem de equacoes estruturais nao vai rodar.
  
  ################# Aqui vamos fazer alteracoes na varia enquanto a matriz
  ################# for NAO definida positiva
  ################# SE teste_covariavianci == F, entao temos que a matriz
  ################ eh NAO DEFINIDA POSITIVA
  
  variaveis_PROBLEMA_MATRIZ<-NULL
  ind<-1
  
  while(teste_covariancia==F){
    
    ## um dos problemas associados com a matriz NAO DEFINIDA positiva eh a 
    ## multicolinearidade. vamos analisar os casos em que ha uma correlacao muito forte
    correl_tirar_var<-cor(dados[,variaveis_do_modelo])
    
    ## colocar zeros na diagonal da matriz
    diag(correl_tirar_var)<-0
    
    ## listar as variaveis com a maior correlacao
    max_correl<-NULL
    correl_tirar_var_aux<-melt(correl_tirar_var)
    max_correl<-correl_tirar_var_aux[order(-correl_tirar_var_aux$value),][1,]
    t2b_var1<-NULL
    t2b_var2<-NULL
    
    ### vamos tirar aquela variavel que tem o menor top2box dai
    t2b_var1<-merge(max_correl, df.top2box, by.x=c("Var1"), by.y=c("variaveis"))
    t2b_var2<-merge(t2b_var1, df.top2box, by.x=c("Var2"), by.y=c("variaveis"))
    
    if(t2b_var2$T2B_result.x[1]==t2b_var2$T2B_result.y[1]){
      
      sorteio_aux<-NULL
      set.seed(1234)
      sorteio_aux<-sample(c(1,2),1)
      t2b_var2$retirar<-as.character(t2b_var2[1,sorteio_aux])
      #################### gravar a variavel acima em um vetor, ########################
      #################### gravar em excel #############################################
      var_retirada_porT2B <- paste(var_retirada_porT2B, "-", as.character(t2b_var2[1,sorteio_aux]))
    }
    else{
      
      t2b_var2$retirar<-ifelse(as.numeric(as.character(t2b_var2$T2B_result.x))<as.numeric(as.character(t2b_var2$T2B_result.y)), 
                               t2b_var2$retirar<-as.character(t2b_var2$Var2),
                               t2b_var2$retirar<-as.character(t2b_var2$Var1))
      
    }
    ## coloca num vetor quais variaveis ficaram de fora
    retirar_var_aux<-NULL
    retirar_var_aux<-t2b_var2$retirar
    variaveis_PROBLEMA_MATRIZ[ind]<-retirar_var_aux
    
    coord_retirar<-NULL
    coord_retirar<-which(variaveis_do_modelo %in% retirar_var_aux)
    variaveis_do_modelo<-variaveis_do_modelo[-coord_retirar]
    
    ## testa se eh positiva definida
    teste_covariancia<-(is.positive.definite(cov(dados[,variaveis_do_modelo])) 
                        & n_obs>length(variaveis_do_modelo))
    
    if(teste_covariancia==F | length(variaveis_PROBLEMA_MATRIZ)>0){
      ind<-ind+1
    }
    else{
      ind<-ind
    }
  }
  
  
  if (length(variaveis_PROBLEMA_MATRIZ) == 0){
    mapeamente_vars_removidas_multicol <- cbind(variaveis= "Nenhuma",motivos="multicolinearidade", momento=2)
  }else{
    mapeamente_vars_removidas_multicol <- cbind(variaveis= variaveis_PROBLEMA_MATRIZ,motivos="multicolinearidade", momento=2)
  }
  
  mapeamente_vars_removidas <- rbind(mapeamente_vars_removidas,mapeamente_vars_removidas_multicol)
  
  ###############################################;
  ## TIRAMOS AS VARIAVEIS QUE COM A MULTICOLINEARIDADE AFETAVAM NA MATRIZ.
  ##  AGORA PRECISAMOS TIRAR AS VARIAVEIS QUE SAIRAM DO MODELO INICIAL
  ## ISSO SO VAI ACONTECER SE FORAM FEITAS ALTERACOES EM FUNCAO
  ## DO TAMANHO DE AMOSTRA OU DE A MATRIZ NAO SER DEFINIDA POSITIVA
  if (ind >1) {
    HS.model2 <- do.remove(HS.model, variaveis_PROBLEMA_MATRIZ)
    
  } else {
    HS.model2 <- HS.model
  }
  
  
  
  
  variaveis_do_modelo_aux_variancia <- NULL
  # Vamos calcular as var?ncias das variaveis que sobraram
  variaveis_do_modelo_aux_variancia <- mapeamente_vars_removidas[,1]
  variaveis_do_modelo_aux_variancia <- setdiff(variaveis_do_modelo, variaveis_do_modelo_aux_variancia)
  #selecionando as linhas da tablea de variancias das variaveis mantidas no moldelo
  variancia_geral_aux <- variancia_geral[variancia_geral[,1] %in% variaveis_do_modelo_aux_variancia,]
  
  variancia_geral_aux <- data.frame(variancia_geral_aux)
  colnames(variancia_geral_aux) <- c("variavel","variancia")
  variancia_geral_aux$concatenar <- NULL
  variancia_geral_aux$concatenar <- paste(paste(paste(variancia_geral_aux$variavel,variancia_geral_aux$variancia,sep="~~"),
                                                variancia_geral_aux$variavel,sep="*"),"  \n")
  #variancia_do_modelo <- paste(variancia_geral_aux$concatenar , collapse = " ")
  #HS.model2 <- paste(HS.model2,variancia_do_modelo)
  
  # vamos armazenar o primeiro modelo gerado do cfa
  cfa_todas_variaveis <- NULL
  cfa_todas_variaveis <- HS.model2
  
  HS.model <- HS.model2
  fit_aux  <- NULL
  fit_aux  <- cfa(model = HS.model,
                  data  = dados,
                  likelihood = "wishart",
                  std.lv = T)
  fit      <- NULL
  modelo_partial <- list()
  modelo_name    <- list()
  ###########################################################################
  #print("inicializando lavaan survey")
  modelo_partial[[1]]<- lavaan.survey(fit_aux, 
                                      svy.df) 
  # print("finalizando lavaan survey - tamanho do arquivo:")
  # print(object.size(modelo_partial))
  ###########################################################################
  razao_chisq<-c()
  aic<-c()
  rmsea<-c()
  cfi<-c()
  ave_min<-c()
  alfa_cronbach<-c()
  valid_discr<-c()
  bic<-c()
  fit_modelo_partial <- NULL
  fit_modelo_partial <- fitMeasures(modelo_partial[[1]])
  razao_chisq[1]<-fit_modelo_partial["chisq"]/fit_modelo_partial["df"]
  aic[1]<-fit_modelo_partial["aic"]
  rmsea[1]<-fit_modelo_partial["rmsea"]
  cfi[1]<-fit_modelo_partial["cfi"]
  alfa_cronbach[1]<-reliability(modelo_partial[[1]])[1,ncol(reliability(modelo_partial[[1]]))]
  ave_min[1]<-min((reliability(modelo_partial[[1]])[5,]))
  bic[1]<-fit_modelo_partial["bic"]
  print("BIC:")
  print(bic)
  print(cfi[1])
  # Validade discriminante
  # Testada atraves da inter-correlacao entre os constructos
  # indices acima de 0.9 podem indicar multicolinearidade entre
  # os constructos (BREI, V."Antecedentes e conseq...", Dissertacao de Mestrado, pg 124)
  mat_validDis<-NULL
  mat_validDis<-inspect(modelo_partial[[1]], "cor.lv")
  diag(mat_validDis)<-NA
  valid_discr[1]<-sum(mat_validDis>0.95 , na.rm=T)
  
  ## utilizando WLS - ADF(Asymptotically Distribution-Free)- Nao supoe nenhuma
  ## nenhuma distribuicao a priori. No nosso caso,as variaveis sao ordinais e nao possuem
  ## distribuicao de suas respostas proximas de uma dist normal
  ## (L. G. GIORDANI, pg 17. "Um modelo de Equacoes Estruturais...", TCC)
  
  ## Porem os metodos de maxima verossimilhanca sao os mais utilizados. Produzem estimativas
  ## dos parametros centradas e consistentes, principalmente a medida que 
  ## se aumenta o tamanho de amostra, alem de serem robustos quanto a violcao da suposicao de normalidade
  ## (S.S. PEREIRA, pg 22. "Modelagem de Equacoes Estruturais no So..", TCC)
  
  ## de acordo com o tutorial do pacote lavaan "The value of the test statistic will 
  ## be closer to the value reported by programs like EQS, LISREL or AMOS,
  ## since they all use the ?Wishart? approach when using the maximum likelihood estimator. 
  ## The program Mplus, on the other hand, uses the ?normal? approach to maximum likelihood estimation."
  ##  http://users.ugent.be/~yrosseel/lavaan/lavaanIntroduction.pdf . pg. 27
  ## Alem disso, permite calcular estatisticas como o BIC e AIC
  
  ### Abaixo temos os resultados da analise confirmatoria para o modelo
  ### completo e original:
  summ_cfa_modelo_comp<-NULL
  estimativas<-NULL
  estimativas<-parameterEstimates(modelo_partial[[1]], rsquare=T)[complete.cases(parameterEstimates(modelo_partial[[1]])),]
  
  # verificar se alguma estimativa eh nao significativa, ou menor que zero 
  est_aux<-NULL
  est_aux<-estimativas[(estimativas$pvalue>0.05
                        & estimativas$op=="=~")| (estimativas$est<0 &
                                                    estimativas$op=="=~"),]
  
  est_aux<-est_aux[complete.cases(est_aux$pvalue),]
  
  ### verificar se nenhum r2 eh maior que 1
  r2_aux<-NULL
  r2_aux<-sum(estimativas$est[estimativas$op=="r2"]>1)
  
  ### verificar se nao ha nenhuma variancia menor que zero
  var_aux<-NULL
  var_aux<-sum(estimativas$est[estimativas$op=="~~"]<0)
  
  ### verificar se tem alguma estimativa padronizada menor que 0.25 Std.all (cuja interpretacao eh igual aos loading da anlise fatorial confirmatoria)
  estimativas_pad<-NULL
  estimativas_pad<-parameterEstimates(modelo_partial[[1]],standardized = T)[complete.cases(parameterEstimates(modelo_partial[[1]])),]
  est_pad_aux<-NULL
  est_pad_aux<-sum(estimativas_pad$std.all[estimativas_pad$op=="=~"]<0.25)
  
  ## vamos usar tambem o RMSEA como um dos criterios de selecao de modelo.
  # De acordo com PEREIRA (2013, pg28  Modelagem de Equacoes Estruturais no Soft R),
  # "eh um dos criterios reconhecidos como mais explicativos na modelagem de estruturas
  # de covariancia, levando em conta o erro de aproximacao na populacao. Ateh um valor 0.05 sao
  # considerados bons; at? 0.08, razoaveis.
  
  ## Outro criterio bastante utilizado eh o de chisq/df
  ## Um valor abaixo de 5 jah eh considerado razoavel.
  
  chisq_df<-fit_modelo_partial["chisq"]/fit_modelo_partial["df"]
  
  if( fit_modelo_partial["rmsea.ci.upper"]<0.08 & fit_modelo_partial["srmr"]<0.08 & 
      fit_modelo_partial["cfi"]>0.8 & nrow(est_aux)<1 & r2_aux<1 & var_aux<1 & 
      est_pad_aux<1 & chisq_df<5 & fit_modelo_partial["df"]>0){
    rodar<-0
  }else{
    rodar<-1
  }
  
  ## sao indices de ajuste geral do modelo construidos a partir
  ## dos residuos obtidos pela diferenca entre a matriz de covariancia
  ## observada e a estimada.  
  ## O RMSEA eh um indice de parcimonia que tenta corrigir as falhas da medida qui-quadrado
  ## o indice penaliza o valor da estatistica de qui-quadrado ao subtrair seu 
  ## valor pelos graus de liberdade. quanto menor melhor
  ## No caso de rodar = 1. Isso quer dizer que o modelo necessita de ajustes:
  
  z<-1
  modelo_name<-list()
  print(z)
  
  # modelo_partial[[1]]<-fit
  modelo_name[[z]] <- HS.model
  
  # print("vai entrar no while da CFA")
  
  while(rodar==1){
    
    if(z==1){
      HS.model<-HS.model
    }else{
      HS.model<-modelofinal
    }
    
    ## primeiro passo eh dividir os models latentes
    latent<-NULL
    latent<-unlist(strsplit(HS.model,
                            split=c("\n")))
    
    latent<-grep("=~", latent, value=T)
    
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
    
    if(z==1){
      # vamos guardar 
      valat_originais<-nome_constructo
    }
    
    #print("lista varia")
    
    # a lista 'varia' tem o nome dos constructos entao
    # vamos nos dedicar a tirar ou colocar variaveis 
    # com objetivo de melhorar o ajuste
    
    # tirar as variaveis nao significativas
    print("ta rodando...se acalme....")
    estimativas<-parameterEstimates(modelo_partial[[z]])[complete.cases(parameterEstimates(modelo_partial[[z]])),]
    est_aux<-NULL
    est_aux<-estimativas[(estimativas$pvalue>0.05
                          & estimativas$op=="=~")| (estimativas$est<0 &
                                                      estimativas$op=="=~")  ,]
    est_aux<-est_aux[complete.cases(est_aux$z),]
    
    retirar<-NULL
    retirar<-NULL
    if (is.na(est_aux$lhs[1])==FALSE) {
      retirar<-cbind(est_aux$lhs, est_aux$rhs)
    }
    
    mapeamento_vars_removidas_baixoimpacto <- NULL
    
    
    mapeamento_vars_removidas_nsignificativas <- NULL
    if(!is.null(retirar)) {
      mapeamento_vars_removidas_nsignificativas <- cbind(variaveis = retirar[,2], 
                                                         motivo = "variavel_nao_significativa_ou_coef_negaivo", momento=z)
    }else{
      mapeamento_vars_removidas_nsignificativas <- cbind(variaveis = "Nenhuma", 
                                                         motivo = "variavel_nao_significativa_ou_coef_negaivo", momento=z)
    }
    
    mapeamente_vars_removidas <- rbind(mapeamente_vars_removidas,mapeamento_vars_removidas_nsignificativas)
    
    i <- 1
    ###########################################################################
    if (typeof(retirar)!="NULL") {
      if (nrow(retirar)>0) {
        for(i in 1:nrow(retirar)) {
          v_l_aux<-NULL
          rhs_aux<-NULL
          v_l_aux<-retirar[i,1]
          rhs_aux<-retirar[i,2]
          
          k<-grep(rhs_aux, varia[[v_l_aux]])
          varia[[v_l_aux]][k]<-NA
          varia[[v_l_aux]]<-varia[[v_l_aux]][complete.cases(varia[[v_l_aux]])]
        } 
      }else{
        retirar<-NULL
      }
    }
    
    
    # print("retirar variaveis nao significativas:ok")
    estimativas2<-NULL
    ###########################################################################
    #### Retirar variaveis cuja padronizacao for menor que 0.10
    estimativas2<-parameterEstimates(modelo_partial[[z]], standardized = T)[complete.cases(parameterEstimates(modelo_partial[[z]])),]
    
    est_aux2<-estimativas2[(estimativas2$std.all<0.10
                            & estimativas$op=="=~"),]
    retirar<-NULL
    if (nrow(est_aux2)!=0) {
      if (is.na(est_aux2$lhs)==FALSE) {
        retirar<-cbind(est_aux2$lhs, est_aux2$rhs)
      } else {
        retirar <- NA
      }
    } else {
      retirar <- NA
    }
    
    
    mapeamento_vars_removidas_baixoimpacto <- NULL
    if(!is.na(retirar)) {
      mapeamento_vars_removidas_baixoimpacto <- cbind(variaveis = retirar[,2], 
                                                      motivo = "variavel_de_baixo_impacto", momento=z)
    }else{
      mapeamento_vars_removidas_baixoimpacto <- cbind(variaveis = "Nenhuma", 
                                                      motivo = "variavel_de_baixo_impacto", momento=z)
    }
    
    mapeamente_vars_removidas <- rbind(mapeamente_vars_removidas,mapeamento_vars_removidas_baixoimpacto)
    
    i <- 1
    if(!is.na(retirar)) {
      for(i in 1:nrow(retirar)) {
        v_l_aux<-NULL
        rhs_aux<-NULL
        v_l_aux<-retirar[i,1]
        rhs_aux<-retirar[i,2]
        
        k<-grep(rhs_aux, varia[[v_l_aux]])
        varia[[v_l_aux]][k]<-NA
        varia[[v_l_aux]]<-varia[[v_l_aux]][complete.cases(varia[[v_l_aux]])]
      } 
    }else{
      retirar<-NULL
    }
    
    #print("retirar variaveis cujo pad for menor que 0.10:ok")
    ###########################################################################
    ## alteracoes para incluir variaveis em constructos
    MI<-NULL
    aval_multi2<-NULL
    if(sum(is.na(parameterestimates(modelo_partial[[z]])$se))/nrow(parameterestimates(modelo_partial[[z]]))<0.5){
      MI<-modificationIndices(modelo_partial[[z]], sort. = T)
      #summary(modelo_partial[[z]])
      ### Ainda o Modification Indices pode nos ajudar a identificar
      ### problemas de multicolinearidade. 
      
      ## Pegar os tres casos que evidenciaram melhorias ao adicionar covariancias 
      multi_aux<-head(subset(MI, mi>200 & op=="~~"),3)
      
      ## colocar num data frame
      multi_aux_a<-NULL
      multi_aux_a<-data.frame(multi_aux$lhs,
                              multi_aux$rhs)
      ## usar o valor das estimativas dos coeficientes nas variaveis latentes
      ## para decidir que variavel tirar. A que tiver menor valor, sai
      estimativas_aux2<-NULL
      maiores_est<-NULL
      estimativas_aux2<-estimativas[estimativas$op=="=~",]
      maiores_est<-data.frame(estimativas_aux2$rhs, estimativas_aux2$est)
      
      if(nrow(maiores_est)>0){
        aval_multi1<-merge(multi_aux_a, maiores_est, by.x="multi_aux.lhs",
                           by.y="estimativas_aux2.rhs",all=F)
      }
      
      if(nrow(aval_multi1)>0 & ncol(aval_multi1)>2){
        
        aval_multi2<-merge(aval_multi1, maiores_est, by.x="multi_aux.rhs",
                           by.y="estimativas_aux2.rhs",all=F)
        
        ## Definir qual variavel tirar:
        if(nrow(aval_multi2)>0){
          aval_multi2$retirar<-ifelse(aval_multi2$estimativas_aux2.est.x>aval_multi2$estimativas_aux2.est.y,
                                      as.character(aval_multi2$multi_aux.rhs),
                                      as.character(aval_multi2$multi_aux.lhs))
          
          for(i in 1:length(aval_multi2$retirar)){
            
            procurar<-NULL
            procurar<-aval_multi2$retirar[i]
            resp<-NULL
            aux<-NULL
            k<-NULL
            
            res <- lapply(varia, function(ch) grep(procurar, ch))
            aux<-sapply(res, function(x) length(x) > 0)
            
            if(sum(aux)>0){
              var_aux<-names(which(aux==T))
              k<-grep(procurar, varia[[var_aux]])
              varia[[var_aux]][k]<-NA
              varia[[var_aux]]<-varia[[var_aux]][complete.cases(varia[[var_aux]])]
            }
          }
        }
      }
    }
    
    mapeamento_vars_removidas_mulitcolMI <- NULL
    if(length(aval_multi2$retirar)>0) {
      mapeamento_vars_removidas_mulitcolMI <- cbind(variaveis = aval_multi2$retirar, 
                                                    motivo = "variavel_de_multicolinearidae_MI",  momento=z)
    }else{
      mapeamento_vars_removidas_mulitcolMI <- cbind(variaveis = "Nenhuma", 
                                                    motivo = "variavel_de_multicolinearidae_MI", momento=z)
    }
    
    mapeamente_vars_removidas <- rbind(mapeamente_vars_removidas,mapeamento_vars_removidas_mulitcolMI)
    
    #print("Tirar multicolinearidade:ok")
    ###########################################################################
    ### Temos que ver ainda se alguma variavel latente nao ficou vazia:
    
    varia<-Filter(length, varia)
    
    var_modelo_revisto <- NULL
    var_modelo_revisto <- data.frame(unlist(varia))
    rownames(var_modelo_revisto) <- NULL
    modelofinal<-NULL
    
    for(i in 1:length(varia)){
      
      va_lt<-NULL
      
      va_lt[i]<-paste(names(varia)[i],"=~")
      
      modelo<-NULL
      var_lat<-(paste(unlist(varia[[i]]),"+"))
      modelo<-paste(var_lat, collapse = " ")
      modelo<-substr(modelo, 1, nchar(modelo)-2)
      modelo<-paste(va_lt[i], modelo,"\n")
      modelofinal<-paste(modelofinal, modelo)
      
    }
    
    #selecionando as linhas da tablea de variancias das variaveis mantidas no moldelo
    variancia_geral_aux_revisto <- NULL
    variancia_geral_aux_revisto <- variancia_geral[variancia_geral[,1] %in% trimws(var_modelo_revisto[,1]),]
    
    variancia_geral_aux_revisto <- data.frame(variancia_geral_aux_revisto)
    colnames(variancia_geral_aux_revisto) <- c("variavel","variancia")
    variancia_geral_aux_revisto$concatenar <- NULL
    variancia_geral_aux_revisto$concatenar <- paste(paste(paste(variancia_geral_aux_revisto$variavel,variancia_geral_aux_revisto$variancia,sep="~~"),
                                                          variancia_geral_aux_revisto$variavel,sep="*"),"  \n")
    #variancia_do_modelo <- paste(variancia_geral_aux_revisto$concatenar , collapse = " ")
    #modelofinal <- paste(modelofinal,variancia_do_modelo)
    
    cfa_GSP_aux<-NULL
    cfa_GSP_aux<-cfa(model = modelofinal,
                     likelihood="wishart",
                     data=dados)
    
    tempo_ponderar_inicio<-NULL
    tempo_ponderar_inicio<-Sys.time()
    
    z<-z+1  
    modelo_partial[[z]]<-lavaan.survey(cfa_GSP_aux, 
                                       svy.df) 
    # print("finalizando lavaan survey - tamanho do arquivo:")
    # print(object.size(modelo_partial))
    ######################################
    print("BIC:")
    bic[z]<-fitMeasures(modelo_partial[[z]], "bic")
    print(bic[z])
    print("CFI:")
    print(fitMeasures(modelo_partial[[z]], "cfi"))
    modelo_name[[z]]<-modelofinal
    estimativas<-NULL
    estimativas<-parameterEstimates(modelo_partial[[z]], rsquare=T)[complete.cases(parameterEstimates(modelo_partial[[z]])),]
    # verificar se alguma estimativa eh nao significativa, ou menor que zero 
    est_aux<-NULL
    est_aux<-estimativas[(estimativas$pvalue>0.05
                          & estimativas$op=="=~")| (estimativas$est<0 &
                                                      estimativas$op=="=~"),]
    est_aux<-est_aux[complete.cases(est_aux$pvalue),]
    ### verificar se nenhum r2 eh maior que 1
    
    r2_aux<-NULL
    r2_aux<-sum(estimativas$est[estimativas$op=="r2"]>1,na.rm=T)
    
    ### verificar se nao ha nenhuma variancia menor que zero
    var_aux<-NULL
    var_aux<-sum(estimativas$est[estimativas$op=="~~"]<0,na.rm=T)
    
    ### verificar se tem alguma estimativa padronizada menor que 0.1
    estimativas_pad<-NULL
    est_pad_aux<-NULL
    chisq_df<-NULL
    estimativas_pad<-parameterEstimates(modelo_partial[[z]],standardized = T)[complete.cases(parameterEstimates(modelo_partial[[z]])),]
    est_pad_aux<-sum(estimativas_pad$std.all[estimativas_pad$op=="=~"]<0.10, na.rm=T)
    fit_modelo_partial <- NULL
    fit_modelo_partial <- fitMeasures(modelo_partial[[z]])
    chisq_df<-fit_modelo_partial["chisq"]/fit_modelo_partial["df"]
    fit<-NULL
    parar<-0
    if(z>2){
      parar<-sum(bic[z-1]==bic[z-2])
    }else{
      parar<-0
    }
    
    if( (fit_modelo_partial["rmsea.ci.upper"]<0.08 & fit_modelo_partial["df"]>0 & nrow(est_aux)<1 & r2_aux<1 & var_aux<1 & est_pad_aux<1 & chisq_df<5) |  (z>30 | parar==1) ){
      rodar<-0
    }else{
      rodar<-1
    }
    adc_aux<-NULL
    aval_multi1<-NULL
    aval_multi2<-NULL
    print(z)
  }
  
  #print("saiu do while da CFA")
  ###########################################################################
  
  if(z>1){
    
    for(i in 1:(z)){
      
      cfi[i]<-fitmeasures(modelo_partial[[i]], "cfi")
      
      
      
    }
    print("*-------------------------------------------------*")
    print("CFI dos modelos de Analise Fatorial Confirmatoria:")
    print(cfi)
    #   ###pegar o modelo com maior cfi
    melhor<-which.max(cfi)
    print("O melhor modelo escolhido foi o de numero:")
    print(melhor)
    print("*-------------------------------------------------*")
    #   ## precisamos tirar os constructos que ficaram de fora:
    #   ## constructos que sobraram:
  }else{
    
    melhor<-1
    
  }
  #print(modelo_partial[[melhor]])
  ## funcao que identifica diferenca entre dois vetores
  outersect <- function(x, y) {
    
    sort(c(setdiff(x, y),
           setdiff(y, x)))
    
  }
  # quais constructos ficaram de fora
  
  diferentes<-outersect(nome_constructo,valat_originais)
  
  ## soh para o caso de detectar mais de 0 constructos que ficaram de fora
  
  if(length(diferentes)>0){
    ## tirando o constructo
    tirar_val<-strsplit(soh_val, split=diferentes, fixed=F)
    
    verificar<-substring(tirar_val[[1]][2], 1, 1)
    
    ### se o primeiro character da segunda parte for '+',
    ## entao precisamos retirar
    if (!is.na(verificar)) {
      if(verificar=="+"){
        
        tirar_val[[1]][2]<-sub('.', '', tirar_val[[1]][2])
        
      }
    }
    
    
    soh_val<-paste(tirar_val[[1]],collapse="")
    ### o constructo retirado pode ainda ser variavel latente de outro constructo;
    ###deve, portanto, na?o estar incluido em soh_val. a variavel diferentes possui
    ###quem deve ser retirado de soh_val.
  }
  
  #constructos <- 
  modelo_completo2<-paste(modelo_name[[melhor]], soh_val, sep = "#")
  
  rodar_modelo2_aux<-NULL
  rodar_modelo2_aux<-sem(modelo_completo2,  
                         likelihood="wishart",
                         data=dados)
  
  survey.fit<-NULL
  survey.fit <- lavaan.survey(rodar_modelo2_aux, 
                              svy.df) 
  # print("finalizando lavaan survey - tamanho do arquivo (survey.fit):")
  # print(object.size(survey.fit))
  while ((survey.fit@Fit@converged)==FALSE) {
    print("n?o convergiu")
    regre_pad_aux<-NULL
    regre_pad_aux<-parameterEstimates(survey.fit,
                                      standardized = T)
    var_neg<-regre_pad_aux[which(regre_pad_aux$op=="~~" & regre_pad_aux$est<0),]
    if (nrow(var_neg)>0) {
      modelo_completo2 <- do.remove(modelo_completo2, var_neg$lhs)
      rodar_modelo2_aux<-NULL
      rodar_modelo2_aux<-sem(modelo_completo2,  
                             likelihood="wishart",
                             data=dados)
      
      survey.fit<-NULL
      survey.fit <- lavaan.survey(rodar_modelo2_aux, 
                                  svy.df) 
      # print("finalizando lavaan survey - tamanho do arquivo (survey.fit):")
      # print(object.size(survey.fit))
    } 
    print((survey.fit@Fit@converged))
    
  }
  ########################################################################################################
  ##### VAMOS INVESTIGAR O MODELO ESTRUTURAL AGORA       #################################################
  ###     AS REGRESSOES                                                                ###################
  ##### Temos que tirar aqueles constructos que deram nao signficativos   ################################
  ##### ou que deram negativos
  
  ## vamos pegar as estimacoes das regressoes
  regre_pad_aux<-NULL
  regre_pad_aux<-parameterEstimates(survey.fit,
                                    standardized = T)[complete.cases(parameterEstimates(survey.fit)),]
  
  ## retornar as var depend e os respectivos constructos nao significatios
  constru_nao_sig<-NULL
  constru_nao_sig<-regre_pad_aux[regre_pad_aux$op=="~" & (regre_pad_aux$pvalue>0.05 | regre_pad_aux$est<0),]
  constru_nao_sig_aux<-cbind(constru_nao_sig$lhs, constru_nao_sig$rhs, constru_nao_sig$est, constru_nao_sig$pvalue)
  constru_nao_sig<-cbind(constru_nao_sig$lhs, constru_nao_sig$rhs)
  
  
  ##objeto auxiliar para pegar as variaveis latentens nao significativas
  constru_nao_sig_aux_reg<-NULL
  constru_nao_sig_aux_reg<-regre_pad_aux[regre_pad_aux$op=="~" & (regre_pad_aux$pvalue>0.05 ),]
  
  ##objeto auxiliar para pegar as variaveis latentens com coeficientes negativos
  constru_nao_negativa_aux_reg<-NULL
  constru_nao_negativa_aux_reg<-regre_pad_aux[regre_pad_aux$op=="~" & (regre_pad_aux$est<0 ),]
  
  # Passar as variaveis latentes nao significativas para lista das variaveis que cairam fora
  
  mapeamento_vars_removidas_regressao_nao_sig <- NULL
  
  if(nrow(constru_nao_sig_aux_reg)>0) {
    mapeamento_vars_removidas_regressao_nao_sig <- cbind(variaveis = constru_nao_sig_aux_reg[,3], 
                                                         motivo = "caiu_na_regressao_nao_significativo",  momento=z)
  }else{
    mapeamento_vars_removidas_regressao_nao_sig <- cbind(variaveis = "Nenhuma", 
                                                         motivo = "caiu_na_regressao_nao_significativo", momento=z)
  }
  
  mapeamente_vars_removidas <- rbind(mapeamente_vars_removidas,mapeamento_vars_removidas_regressao_nao_sig)
  
  # Passar as variaveis latentes negativas para lista das variaveis que cairam fora
  
  mapeamento_vars_removidas_regressao_negativas <- NULL
  
  if(nrow(constru_nao_negativa_aux_reg)>0) {
    mapeamento_vars_removidas_regressao_negativas <- cbind(variaveis = constru_nao_negativa_aux_reg[,3], 
                                                           motivo = "caiu_na_regressao_negativo", momento=z)
  }else{
    mapeamento_vars_removidas_regressao_negativas <- cbind(variaveis = "Nenhuma", 
                                                           motivo = "caiu_na_regressao_negativo", momento=z)
  }
  
  mapeamente_vars_removidas <- rbind(mapeamente_vars_removidas,mapeamento_vars_removidas_regressao_negativas)
  
  ### verifcar se precisa algum ajuste no modelo estrutural
  ver_ajuste<-nrow(constru_nao_sig)>0
  contador<-1
  qtos_prob<-NULL
  
  while(ver_ajuste==T){
    
    if(contador==1){
      
      modelo_reg<-soh_val
      
    }else{
      
      modelo_reg<-modelo_reg
    }
    
    pegar_reg<-NULL
    pegar_reg<-unlist(strsplit(modelo_reg,
                               split=c("\n")))
    pegar_reg2<-NULL
    pegar_reg2<-grep("~", pegar_reg, value=T)
    
    pegar_reg3<-list()
    pegar_reg4<-c()
    
    for(i in 1:length(pegar_reg2)){
      pegar_reg3[i]<-strsplit(pegar_reg2[i],
                              split=c("~"))
      pegar_reg4[i]<-trimws(pegar_reg2[[i]][1])
      pegar_reg4<-trimws(pegar_reg4)
    }
    
    #print("primeira etapa da separacao variaveis")
    
    # separar os preditores dentro de cada constructo
    
    pegar_reg5<-list()
    var_dep_reg<-c()
    
    for(i in 1:length(pegar_reg3)){
      
      ## vamos dividir em cada elemento da parte dos preditores 
      
      pegar_reg5[i]<-strsplit(pegar_reg3[[i]][2],
                              split=c("\\+"))
      
      var_dep_reg[i]<-pegar_reg3[[i]][1]  
    }
    
    #print("separar os preditores dentro de cada construco")
    ## lista das regressoes
    names(pegar_reg5)<-trimws(var_dep_reg)
    
    if(contador==1){
      ## vamos pegar as estimacoes das regressoes
      regre_pad_aux<-NULL
      regre_pad_aux<-parameterEstimates(survey.fit,
                                        standardized = T)[complete.cases(parameterEstimates(survey.fit)),]
      
      ## retornar as var depend e os respectivos constructos nao significatios
      constru_nao_sig<-regre_pad_aux[regre_pad_aux$op=="~" & (regre_pad_aux$pvalue>0.05 | regre_pad_aux$est<0),]
      constru_nao_sig<-cbind(constru_nao_sig$lhs, constru_nao_sig$rhs)
    }
    
    
    for(i in 1:nrow(constru_nao_sig)){
      
      ver_aux<-NULL
      rhs_aux<-NULL
      k<-NULL
      ver_aux<-constru_nao_sig[i,1]
      rhs_aux<-constru_nao_sig[i,2]
      
      tab_aux<-NULL
      tab_aux<-lapply(pegar_reg5, function(ch) grep(rhs_aux, ch))
      
      # verificar se a variavel nao esta em nenhum constructo
      k<-grep(rhs_aux, pegar_reg5[[ver_aux]])
      pegar_reg5[[ver_aux]][k]<-NA
      pegar_reg5[[ver_aux]]<-pegar_reg5[[ver_aux]][complete.cases(pegar_reg5[[ver_aux]])]
    }
    
    modelo_reg<-NULL
    
    for(i in 1:length(pegar_reg5)){
      #i <- 2
      va_lt<-NULL
      
      va_lt[i]<-paste(names(pegar_reg5)[i],"~")
      
      modelo<-NULL
      var_lat<-(paste(unlist(pegar_reg5[[i]]),"+"))
      modelo<-paste(var_lat, collapse = " ")
      modelo<-substr(modelo, 1, nchar(modelo)-2)
      modelo<-paste(va_lt[i], modelo,"\n")
      modelo_reg<-paste(modelo_reg, modelo)
    }
    
    
    ############# vamos tirar aquelas regressoes cujas variaveis independentes
    ############# ficaram vazias
    
    
    modelo_reg_aux<-strsplit(modelo_reg,
                             split=c("\\n"))
    
    for(i in 1:length(modelo_reg_aux[[1]])){
      modelo_reg_aux2 <- NULL
      modelo_reg_aux2<-strsplit(trimws(modelo_reg_aux[[1]][i]),
                                split=c("\\~"))
      
      
      if(length(modelo_reg_aux2[[1]])<2){
        modelo_reg_aux[[1]][i]<-NA
      }
      
    }
    modelo_reg_aux <- modelo_reg_aux[[1]][complete.cases(modelo_reg_aux)]
    modelo2<-NULL
    modelo_reg2 <- NULL
    for (i in 1:length(modelo_reg_aux)) {
      var_lat_aux2 <- NULL
      var_lat_aux2<-(paste(unlist(modelo_reg_aux[i]),"\n"))
      modelo2<-paste(var_lat_aux2, collapse = " ")
      modelo_reg2<-paste(modelo_reg2, modelo2)
    }
    
    
    modelo_reg <- modelo_reg2
    
    modelo_completo3<-NULL
    modelo_completo3<-paste(modelo_name[[melhor]], modelo_reg)
    
    rodar_modelo3_aux<-sem(modelo_completo3,  
                           likelihood="wishart",std.lv=T,
                           data=dados)
    
    tempo_ponderar_inicio<-NULL
    tempo_ponderar_inicio<-Sys.time()
    
    survey.fit2 <- lavaan.survey(rodar_modelo3_aux, 
                                 svy.df)
    #print("finalizando lavaan survey - tamanho do arquivo (survey.fit2):")
    # print(object.size(survey.fit2))
    ######################################
    ## vamos pegar as estimacoes das regressoes
    regre_pad_aux<-NULL
    regre_pad_aux<-parameterEstimates(survey.fit2,
                                      standardized = T)[complete.cases(parameterEstimates(survey.fit)),]
    
    ## retornar as var depend e os respectivos constructos nao significatios
    constru_nao_sig<-NULL
    constru_nao_sig<-regre_pad_aux[regre_pad_aux$op=="~" & (regre_pad_aux$pvalue>0.05 | regre_pad_aux$est<0),]
    constru_nao_sig<-constru_nao_sig[complete.cases(constru_nao_sig),]  
    constru_nao_sig<-(cbind(constru_nao_sig$lhs, constru_nao_sig$rhs))
    
    ver_ajuste<-(nrow(constru_nao_sig)>0 & contador<2)
    qtos_prob[contador]<-nrow(constru_nao_sig)
    print(qtos_prob[contador])
    contador<-contador+1
    #print(contador)
    
  }
  
  ### Modelo final
  var_fora <- NULL
  if(contador>1){
    rodar_modelo3_aux<-sem(modelo_completo3,  
                           likelihood="wishart",
                           data=dados)
    var_fora <- modelo_completo3
    modelo_final_final<- lavaan.survey(rodar_modelo3_aux, 
                                       svy.df) 
    # print("finalizando lavaan survey - tamanho do arquivo (modelo_final_final):")
    # print(object.size(modelo_final_final))
  }else{
    rodar_modelo2_aux<-sem(modelo_completo2,  
                           likelihood="wishart",
                           data=dados)
    var_fora <- modelo_completo2
    modelo_final_final<- lavaan.survey(rodar_modelo2_aux, 
                                       svy.df) 
    # print("finalizando lavaan survey - tamanho do arquivo (modelo_final_final):")
    # print(object.size(modelo_final_final))
    
  }
  
  banco_calc_indireto<-NULL
  banco_calc_indireto<-parameterEstimates(modelo_final_final, standardized = T)[parameterEstimates(modelo_final_final,standardized = T)$op == "~" & parameterEstimates(modelo_final_final,standardized = T)$std.all > 0.1,]
  
  if(!(eh_CES) && var_CB %in% banco_calc_indireto$rhs){
    
    
    estimativa_CB<-NULL
    estimativa_CB<-banco_calc_indireto[banco_calc_indireto$rhs==var_CB,"std.all"]
    
    
    # vamos descobrir quais variaveis estao tanto na regressao para Custo beneficio
    # quanto estao na de satisfacao geral
    
    vetor_var_cb<-NULL
    vetor_var_cb<-banco_calc_indireto[banco_calc_indireto$lhs==var_CB,"rhs"]
    vetor_var_sat<-NULL
    vetor_var_sat<-banco_calc_indireto[banco_calc_indireto$lhs==var_SAT,"rhs"]
    
    variaives_em_comum<-NULL
    variaives_em_comum<-intersect(vetor_var_cb,vetor_var_sat)
    if (length(variaives_em_comum)!=0) {
      if (length(vetor_var_cb)!=0) {
        
        
        banco_efeitos_totais<-NULL
        banco_efeitos_totais<-banco_calc_indireto[banco_calc_indireto$lhs==var_CB,c("rhs","std.all")]
        banco_efeitos_totais<-banco_efeitos_totais[banco_efeitos_totais$rhs %in% variaives_em_comum,]
        banco_efeitos_sat <- banco_calc_indireto[banco_calc_indireto$lhs==var_SAT & banco_calc_indireto$rhs %in% variaives_em_comum,]
        banco_efeitos_totais$efeito_CB<-estimativa_CB
        banco_efeitos_totais$efeitos_indiretos<-banco_efeitos_totais$std.all*banco_efeitos_totais$efeito_CB
        banco_efeitos_totais$total<- banco_efeitos_sat[,"std.all"]+banco_efeitos_totais$efeitos_indiretos
        var_somente_custobenificio<- NULL
        var_somente_custobenificio<- vetor_var_cb[vetor_var_cb %in% vetor_var_sat==F] 
        
        #### se o tamanho das variaveis em comum for igual ao das somente em custo beneficio, nao temos variaveis apenas em CB:
        if (length(var_somente_custobenificio)!=0) {
          
          #### variaveis que estavam no Custo beneficico mas nao estavam no de Satisfacao Geral
          
          
          ### calcular o efeito dessa variavel em funcao do custo beneficio
          
          banco_calc_Sat_via_CB<-NULL
          banco_calc_Sat_via_CB<-banco_calc_indireto[banco_calc_indireto$rhs %in% var_somente_custobenificio & banco_calc_indireto$lhs==var_CB,c("rhs","std.all")]
          
          banco_calc_Sat_via_CB$efeito_CB<-estimativa_CB
          banco_calc_Sat_via_CB$efeitos_indiretos<-estimativa_CB*banco_calc_Sat_via_CB$std.all
          banco_calc_Sat_via_CB$total<-banco_calc_Sat_via_CB$efeitos_indiretos
          
          banco_efeitos_totais<-rbind(banco_efeitos_totais,banco_calc_Sat_via_CB)
          
        }
      }
    }
    
    
    
    
    
  } else {
    #######################################################
    ####  A ideia agora eh estimar os scores de cada individuo para cada variavel latentne
    
    banco_para_cargas_original<-NULL
    banco_para_cargas<-NULL
    banco_guarda_scores<-NULL
    banco_cargas_segunda_ordem<-NULL
    banco_para_cargas_original<-parameterestimates(modelo_final_final)
    banco_para_cargas<-banco_para_cargas_original[banco_para_cargas_original$op=="=~",]
    
    ### Num primeiro momento, vamos calcular apenas aquelas
    ### variaveis de primeira ordem.
    banco_cargas_segunda_ordem<-banco_para_cargas[banco_para_cargas$lhs == segunda_ordem,]
    banco_para_cargas<-banco_para_cargas[banco_para_cargas$lhs != segunda_ordem,]
    
    # quais variaveis latentes que sobraram
    var_latentes_scores<-NULL
    var_latentes_scores<-unique(banco_para_cargas$lhs)
    cont_cargas<-NULL
    cont_cargas<-1
    
    for(cont_cargas in 1:length(var_latentes_scores)){
      ## pegando a "cont_cargas" variavel latente
      var_latente_aux_score<-NULL
      var_latente_aux_score<-var_latentes_scores[cont_cargas]
      
      banco_para_cargas_aux<-NULL
      banco_para_cargas_aux<-banco_para_cargas[banco_para_cargas$lhs %in% var_latente_aux_score, c("rhs","est")]
      banco_para_cargas_aux<-banco_para_cargas_aux[order(banco_para_cargas_aux$rhs),]
      
      ## Um banco de dados soh para cargas
      data_aux_cargas<-NULL
      data_aux_cargas<-dados
      
      ### Multiplicacao das cargas pelos valores das variaveis observadas
      
      mult_cargas<-NULL
      #n_ind_calc_carga<-NULL
      #n_ind_calc_carga<-nrow(data_aux_cargas)
      
      #banco_ordenado<-NULL
      #banco_ordenado<-with(data_aux_cargas,  banco_para_cargas_aux[order(rhs) , ])
      
      
      
      mult_cargas<-t(t(matrix(banco_para_cargas_aux[,'est'])) %*% t(as.matrix(data_aux_cargas[,sort(banco_para_cargas_aux$rhs)])))
      banco_guarda_scores_aux<-NULL
      banco_guarda_scores_aux<-mult_cargas/sum(banco_para_cargas_aux[,2])
      
      banco_guarda_scores<-cbind(banco_guarda_scores,banco_guarda_scores_aux)
      colnames(banco_guarda_scores)<-var_latentes_scores
      #with(df,  df[order(sortbythiscolumn) , ])
      
      ## Calculando os escores de segunda ordem
      
      
      
    }
    banco_cargas_segunda_ordem<-with(banco_cargas_segunda_ordem,  banco_cargas_segunda_ordem[order(rhs) , ])
    
    segunda_ordem_SCORE<-NULL
    segunda_ordem_SCORE<-t(t(matrix(banco_cargas_segunda_ordem[,"est"])) %*% t(as.matrix(banco_guarda_scores[,order(colnames(banco_guarda_scores))])))
    banco_guarda_scores<-cbind(trimws(as.character(dados$ResponseId)), banco_guarda_scores,SCORE_FINAL=segunda_ordem_SCORE/sum(banco_cargas_segunda_ordem[,"est"]))
    
    colnames(banco_guarda_scores)<-c("ResponseId",var_latentes_scores, "Score_final")
  }
  
  
  
  
  
  
  #### Alem disso, vamos rodar o modelo completasso
  modelo_origem_pra_finalera<-NULL
  modelo_origem_pra_finalera<-paste(cfa_todas_variaveis, soh_val)
  SEM_modelo_origem<-NULL
  SEM_modelo_origem<-sem(modelo_origem_pra_finalera,  
                         likelihood="wishart",std.lv=T,
                         data=dados)
  SEM_modelo_origem_final<-NULL
  SEM_modelo_origem_final<- lavaan.survey(SEM_modelo_origem, 
                                          svy.df)
  
  #print("finalizando lavaan survey - tamanho do arquivo (SEM_modelo_origem_final):")
  #print(object.size(SEM_modelo_origem_final))
  
  ###########################################################################
  ####                      EXPORTAcaO DO RESULTADO                      ####
  ###########################################################################
  ## Definir diretorio para onde serao exportados os resultados
  #setwd(export)
  
  ## As estimativas dos parametros do modelo final:
  est_modelo_final<-NULL
  est_modelo_final<-parameterEstimates(modelo_final_final, 
                                       rsquare = T, standardized = T)
  
  # criar uma variavel id para ordenara a tabela depois
  est_modelo_final$ID<-seq(1,nrow(est_modelo_final))
  ## de para para pegar o label das variaveis
  DePara_var<-merge(est_modelo_final, var_labels, by.x="rhs", by.y="codigo", all.x=T, sort = F)
  names(DePara_var)[names(DePara_var) == 'label'] <- 'label_rhs'
  DePara_var2<-merge(DePara_var, var_labels, by.x="lhs", by.y="codigo", all.x=T, sort = F)
  names(DePara_var2)[names(DePara_var2) == 'label'] <- 'label_lhs'
  ## Ordernar de acordo com ID e ordenar as colunas
  DePara_var2<-DePara_var2[order(DePara_var2$ID),]
  est_modelo_final2<-NULL
  est_modelo_final2<-DePara_var2[order(DePara_var2$op),]
  est_modelo_final2<-est_modelo_final2[,c("lhs","label_lhs","op","rhs","label_rhs","est","se","z","pvalue","ci.lower","ci.upper","std.lv","std.all","std.nox","ID")]
  est_modelo_final2<-est_modelo_final2[order(est_modelo_final2$ID),]
  
  est_modelo_final2$label_rhs<-as.character(est_modelo_final2$label_rhs)
  est_modelo_final2$label_lhs<-as.character(est_modelo_final2$label_lhs)
  
  est_modelo_final2$label_rhs <- ifelse(is.na(est_modelo_final2$label_rhs),
                                        est_modelo_final2$rhs, 
                                        est_modelo_final2$label_rhs)
  
  est_modelo_final2$label_lhs <- ifelse(is.na(est_modelo_final2$label_lhs), 
                                        est_modelo_final2$lhs, 
                                        est_modelo_final2$label_lhs)
  
  est_modelo_final2<-est_modelo_final2[(est_modelo_final2$op=="~" | est_modelo_final2$op=="=~" | est_modelo_final2$op=="r2"),]
  ## colocar colunas com as interpretacoes dos resultados
  ## Interpretacoees baseadas en https://www.lume.ufrgs.br/bitstream/handle/10183/133719/000986124.pdf?sequence=1
  ## paginas 51 e posteriores:
  
  ## Outra fonte muito boa de interpretacao de resultados:
  
  
  # "Factor loadings can be interpreted like a regression coefficient. 
  # For example, the first parameter says that for each unit
  # increase in the latent visual ability 
  # (since we standardized latent factors, this means for each 1SD increase),
  # the model predicts a .90-unit increase in x1. Because we included
  # standardized=TRUE in the command, standardized parameter 
  # estimates as well; the column called "std.all" is the one people 
  # typically think of as "standardized regression coefficients" 
  # (often reported as BETA in regression software). 
  # The full output actually includes two additional 
  # columns of standardized coefficients,
  # but I omitted them here to save room and because they're rarely reported."
  
  #http://www.understandingdata.net/2017/03/22/cfa-in-lavaan/
  # Rose Hartman, PhD
  
  
  est_modelo_final_val_aux<-est_modelo_final2[(est_modelo_final2$op=="=~"),]
  est_modelo_final_reg_aux<-est_modelo_final2[(est_modelo_final2$op=="~"),]
  est_modelo_final_reg_val<-rbind(est_modelo_final_val_aux,
                                  est_modelo_final_reg_aux[which(est_modelo_final_reg_aux$std.all>0.1),])
  est_modelo_final_reg_val<-merge(est_modelo_final_reg_val,
                                  df.top2box,
                                  by.x="rhs",
                                  by.y="variaveis", all.x=T)
  
  est_modelo_final_reg_val<-est_modelo_final_reg_val
  
  r2_bancofim<-est_modelo_final2[est_modelo_final2$op=="r2",]
  r2_bancofim<-r2_bancofim[,colSums(is.na(r2_bancofim))<nrow(r2_bancofim)]
  r2_bancofim_v1<-r2_bancofim[!duplicated(as.list(r2_bancofim))]
  
  est_modelo_final3<-est_modelo_final2[est_modelo_final2$op!="r2",]
  
  ###########################################################################
  # impacto
  
  est_modelo_impacto_aux<-est_modelo_final2[(est_modelo_final2$op=="~" | est_modelo_final2$op=="=~") & est_modelo_final2$std.all>0.1,]
  
  est_modelo_impacto<-merge(est_modelo_impacto_aux,
                            df.top2box,
                            by.x="rhs",
                            by.y="variaveis", all.x=T)
  
  #est_modelo_impacto<-est_modelo_impacto[c("lhs", "label_lhs", "op", "rhs", "label_rhs","std.all","T2B_result")]
  est_modelo_impacto$T2B_result<-as.numeric(as.character(est_modelo_impacto$T2B_result))
  
  ################################################################################################################################################################################################################
  ### Exportar para um excel
  ################################################################################################################################################################################################################
  
  # Create an Excel workbook. 
  #setwd(export)
  # vamos criar uma nova pasta das quebras
  data <- date()
  hour <- (strsplit(data," ")[[1]])[4]%>% 
    stringr::str_replace_all(.,":","")
  hour <- as.numeric(hour)+sample(1:999, 1)
  fileName <- paste(Estudo,
                    paste(nivel_label,
                          paste(indice,
                                paste(hour,"comN500-sem_outliers.xlsx", sep=""),
                                sep="_"),
                          sep="_"),
                    sep="_")  %>%
    stringr::str_replace_all(.,"/","")  %>%
    stringr::str_replace_all(.," ","")  %>%
    stringr::str_replace_all(.,"-","")  
  fileXls <- fileName
  unlink(fileXls, recursive = FALSE, force = FALSE)
  exc <- loadWorkbook(fileXls, create = TRUE)
  
  #---- Estimativas na planilha estimativas ----
  createSheet(exc,'Cargas')
  saveWorkbook(exc)
  input <- est_modelo_final_reg_val
  writeWorksheet(exc, input, sheet ='Cargas', startRow = 1, startCol = 2)
  saveWorkbook(exc)
  
  #---- Estatisticas resumo modelo ----
  resumo_modelo<-fitmeasures(modelo_final_final, c("rmsea", "cfi", "df"))
  createSheet(exc,'EstatResumo')
  saveWorkbook(exc)
  input <- as.data.frame(cbind(round(resumo_modelo,4), 
                               rownames(resumo_modelo)))
  input<-rbind(input,n=round(n_obs,1), z=melhor)
  
  input$estat<-rownames(input)
  writeWorksheet(exc, input, sheet ='EstatResumo',header=T,rownames = T, startRow = 1, startCol = 2)
  saveWorkbook(exc)
  #---- Modelo_FINAL_medidas resum ----
  createSheet(exc,'Mod_FINAL_FITMEASUR')
  input<-NULL
  input_aux<-NULL
  saveWorkbook(exc)
  input_aux <- fitMeasures(modelo_final_final)
  input_aux<-data.frame(input_aux)
  input<-data.frame(medidas=rownames(input_aux),valores=format(input_aux,scientific=FALSE))
  
  
  writeWorksheet(exc, input,  sheet ='Mod_FINAL_FITMEASUR',header=T,rownames = T, startRow = 1, startCol = 2)
  saveWorkbook(exc)
  
  #---- R squared ----
  createSheet(exc,'Rsquared')
  saveWorkbook(exc)
  input <- r2_bancofim_v1
  writeWorksheet(exc, input, sheet ='Rsquared',header=T,rownames = T, startRow = 1, startCol = 2)
  saveWorkbook(exc)
  
  #---- Reliability CFA ----
  createSheet(exc,'ReliabilityCFA')
  saveWorkbook(exc)
  rel_aux<-reliability(modelo_partial[[z]])
  input <- data.frame(Estat=row.names(rel_aux),rel_aux) 
  writeWorksheet(exc, input, 
                 sheet ='ReliabilityCFA',header=T,rownames = T, startRow = 1, startCol = 2)
  
  if(!(eh_CES)){
    for(i in 1:length(segunda_ordem)){
      input2<-NULL
      input3<-NULL
      input2<-data.frame(reliabilityL2(modelo_final_final,segunda_ordem[i]))
      input3<-data.frame(cbind(CR=row.names(input2),segunda_ordem[i],est=(input2[,1])))
      input3$est<-as.numeric(as.character(input3$est))
      
      
      writeWorksheet(exc, input3, 
                     sheet ='ReliabilityCFA',header=T,rownames = T, startRow = nrow(input)+
                       (7*i), startCol = 2)
    }
    
    modelo_partial = NULL
    saveWorkbook(exc)
  }
  #---- Impactos ----
  createSheet(exc,'Impactos')
  saveWorkbook(exc)
  input<-NULL
  input <- est_modelo_impacto[order(est_modelo_impacto$lhs),]
  writeWorksheet(exc, input, sheet ='Impactos',header=T,rownames = T, startRow = 1, startCol = 2)
  saveWorkbook(exc)
  
  #---- Efeitos Indiretos e Totais via Custo Beneficio ----
  if(!(eh_CES)){
    if (length(variaives_em_comum)!=0) {
      if (length(vetor_var_cb)!=0) {
        colnames(banco_efeitos_totais)<-c('variavel','efeito_orig', 'efeito_CB', 'efeito_indireto', 'efeito_final')
        createSheet(exc,'Efeitos_Totais_SATGER')
        saveWorkbook(exc)
        input<-NULL
        
        input <- banco_efeitos_totais
        
        writeWorksheet(exc, input, sheet ='Efeitos_Totais_SATGER',header=T,rownames = T, startRow = 1, startCol = 2)
        saveWorkbook(exc)
      }
    }
  }
  #---- Modelo utilizado ----
  createSheet(exc,'ModeloEscolh')
  saveWorkbook(exc)
  input<-NULL
  input <- modelo_completo2
  writeWorksheet(exc, input, sheet ='ModeloEscolh',header=T,rownames = T, startRow = 1, startCol = 2)
  saveWorkbook(exc)
  #---- VarFora ----
  createSheet(exc,'VarFora')
  saveWorkbook(exc)
  input<-NULL
  mapeamente_vars_removidas<-mapeamente_vars_removidas[- grep("Nenhuma", mapeamente_vars_removidas[,1]),]
  input <- mapeamente_vars_removidas
  writeWorksheet(exc, input, sheet ='VarFora',header=T,rownames = T, startRow = 1, startCol = 2)
  #writeWorksheet(exc, model_completo_aux, sheet ='VarFora',header=T,rownames = T, startRow = 1, startCol = 3)
  saveWorkbook(exc)
  #---- Modelo_ORIGEM ----
  createSheet(exc,'Modelo_ORIGEM')
  input<-NULL
  saveWorkbook(exc)
  input <- parameterestimates(SEM_modelo_origem_final,standardized = T)
  writeWorksheet(exc, input, sheet ='Modelo_ORIGEM',header=T,rownames = T, startRow = 1, startCol = 2)
  saveWorkbook(exc)
  
  #---- Modelo_ORIGEM_medidas resum ----
  createSheet(exc,'Mod_ORIGEM_FITMEASUR')
  input<-NULL
  input_aux<-NULL
  saveWorkbook(exc)
  input_aux <- fitMeasures(SEM_modelo_origem_final)
  input_aux<-data.frame(input_aux)
  input<-data.frame(medidas=rownames(input_aux),valores=format(input_aux,scientific=FALSE))
  
  
  writeWorksheet(exc, input,  sheet ='Mod_ORIGEM_FITMEASUR',header=T,rownames = T, startRow = 1, startCol = 2)
  saveWorkbook(exc)
  if ((eh_CES)) {
    #----------- Salvando Scores por individuo------
    createSheet(exc,'Scores_ind')
    input<-NULL
    input_aux<-NULL
    saveWorkbook(exc)
    input_aux <- banco_guarda_scores
    input_aux<-data.frame(input_aux)
    input<-data.frame(medidas=rownames(input_aux),valores=format(input_aux,scientific=FALSE))
    writeWorksheet(exc, input,  sheet ='Scores_ind',header=T,rownames = T, startRow = 1, startCol = 2)
    saveWorkbook(exc)
    
  }
  
  drop_upload(fileName, path = IES)
  #------------------
  
  
  
}

#save.image("Satisfaction_2018.RData")
#history("Satisfaction_2018.Rhistory")

#https://www.dropbox.com/s/1yhpztfyx4z4q5v/UNIFACS2_Satisfaction_2019.sav?dl=0
#https://www.dropbox.com/s/p1o5favhy86szg4/UNIFACS3_small_Satisfaction_2019.sav?dl=0
IES <<- "UNIFACS2"
n_quebra <- 1
pular <- 1
quebras <- c("Finished")
cod_dropbox <- "1yhpztfyx4z4q5v"
eh_CES<-FALSE
if (eh_CES) {
  model_completo <- '
  apoio=~CES_27+CES_34+CES_35+CES_36+CES_37+CES_38+CES_39+CES_41+CES_43+CES_45
  formacao=~CES_01+CES_02+CES_03+CES_04+CES_07+CES_10+CES_11
  facilities=~CES_46+CES_47+CES_48+CES_49+CES_50+CES_51+CES_52
  first=~CES_15+CES_17+CES_18+CES_19+CES_20+CES_21
  portal=~CES_31+CES_32
  CES=~apoio+formacao+facilities+first+portal
  #
  '
  var_CB<-NULL
  var_SAT<-NULL
  constructs <- c("apoio","formacao","facilities","first","portal")
  var_dep_eq <-NULL
  segunda_ordem <-  c("CES")
  var_lonely    <-  NULL
  
  quebras_excluidas   <- NULL
  var_retirada_porT2B <- NULL
} else {
  
  ########### definicao das specs do modelo (CES) #########
  
  #################### modelo para satisfaction (SEM): ################
  model_completo <- '
  classroom=~INF_02_03+INF_02_04+INF_02_05+INF_02_06+INF_02_07
  it_labs=~INF_03_01+INF_03_02+INF_03_03+INF_03_04+INF_03_06
  specific_labs=~INF_04_01+INF_04_02+INF_04_03+INF_04_06
  library=~INF_05_02+INF_05_04+INF_05_10+INF_05_09+INF_05_07+INF_05_21+INF_05_20
  blackboard=~INF_06_01+INF_06_02+INF_06_04+INF_06_07+INF_06_08
  program=~PROG_02_01+PROG_02_03+PROG_02_04+PROG_02_05
  disc_online=~ELEAR_02_01+ELEAR_02_02+ELEAR_02_03+ELEAR_02_04+ELEAR_02_05+ELEAR_02_06
  faculty=~PROF_02_01+PROF_02_02+PROF_02_03+PROF_02_04+PROF_02_05+PROF_02_06+PROF_02_07+PROF_02_08+PROF_02_15
  coord=~COORD_02_01+COORD_02_02+COORD_02_03+COORD_02_04
  financial_serv=~FSERV_01_01+FSERV_01_02+FSERV_01_03+FSERV_01_04
  std_services=~SSERV_01_01+SSERV_01_02+SSERV_01_03+SSERV_01_04
  comunic=~COMM_01_01+COMM_01_02+COMM_01_03+COMM_01_04+COMM_01_05+COMM_01_06+COMM_01_07
  image=~IMG_01_01+IMG_01_02+IMG_01_03
  campus=~INF_01_01+INF_01_02+INF_01_04+INF_01_06+INF_01_09+INF_01_11+INF_01_13+INF_01_15+INF_01_12+INF_01_08
  call=~CALL_01_01+CALL_01_02+CALL_01_03+CALL_01_04
  employ=~EMPL_02_01+EMPL_02_02+EMPL_02_03+EMPL_02_04+EMPL_02_05+EMPL_02_06
  inter=~INTL_01_01+INTL_01_02+INTL_01_03+INTL_01_04+INTL_01_05+INTL_01_07
  infra=~campus+classroom+it_labs+specific_labs+library+blackboard
  services=~financial_serv+std_services+call+employ+inter
  #
  SAT_00_02~infra+faculty+coord+program+disc_online+comunic+services+image
  SAT_00_01~infra+faculty+coord+program+disc_online+comunic+services+image+SAT_00_02
  NPS_01_01~SAT_00_01+image'
  
  ## definir qual eh codigo da variavel custo beneficio (var_CB)
  var_CB<-"SAT_00_02"
  # 
  
  ## definir qual eh codigo da variavel Satisfacao Geral (var_SAT)
  
  var_SAT<-"SAT_00_01"
  # 
  
  # 
  constructs <- c("classroom","it_labs","specific_labs","library","blackboard","program","disc_online","faculty","coord","financial_serv","std_services","comunic","image","campus","call","employ","inter")
  
  var_dep_eq <-c("SAT_00_02","SAT_00_01","NPS_01_01")
  # 
  segunda_ordem <-  c("infra","services")
  var_lonely    <-  NULL
  
  quebras_excluidas   <- NULL
  var_retirada_porT2B <- NULL
}

for (i in 1:100) {
  do.boot(IES = IES,
          quebras, n_quebra, pular,
          model_completo, constructs, var_dep_eq, segunda_ordem, var_lonely,
          cod_dropbox)
}
#Rprof(tf <- "rprof.log", memory.profiling=TRUE)

#Rprof(NULL)
#summaryRprof(tf)
