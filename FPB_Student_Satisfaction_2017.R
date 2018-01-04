install.packages("lavaan")

install.packages("foreign")

install.packages("semTools")

install.packages("rJava")

install.packages("XLConnect")

install.packages("survey")

install.packages("tidyverse")

install.packages("doParallel")

install.packages("rdrop2")

install.packages("RCurl")

install.packages("httpuv")

install.packages("matrixcalc")

install.packages("reshape2")

library(httpuv)
library(lavaan)
library(foreign)
library(semTools)
library(rJava)
library(XLConnect)
library(survey)
library(tidyverse)
library(doParallel)
library(rdrop2)
library(RCurl)
library(matrixcalc)
library(reshape2)

# Diretorio de onde baixar o banco de dados
setwd("/home/gabrielpehls/sem")

#memory.limit(6384)
########################################################################################
########################------  Funções utilizadas  ------##############################
########################################################################################
############---- Download do spss e dos arquivos para upload no dropbox ----############
dl_from_dropbox <- function(x, key) {
  require(RCurl)
  bin <- getBinaryURL(paste0("https://dl.dropboxusercontent.com/s/", key, "/", x),
                      ssl.verifypeer = FALSE)
  con <- file(x, open = "wb")
  writeBin(bin, con)
  close(con)
  message(noquote(paste(x, "read into", getwd())))
}
dl_from_dropbox("FPB_Student%20Satisfaction%202017_DB_V1_Gradua%C3%A7%C3%A3o.sav","1pkzr0q6g3x3pg9")##################---- Arquivos para autorizar upload no dropbox -----##################
##############-----Devem ser atualizados para outra conta do dropbox------##############
dl_from_dropbox(".httr-oauth","1j0froxvx83quxd")
dl_from_dropbox("token.rds","169we6hudfy83u8")
token <- readRDS("token.rds")
drop_acc(dtoken = token)
#########################---- para fazer upload no dropbox ----#########################
#drop_upload('Law_resultados_SEM.xlsx')
########################################################################################
difference <- list()
getDifference <- function (difference) {
  formulas <- unlist(strsplit(modelo_name[[melhor]],
                              split=c("\n")))
  parametros <- NULL 
  for (i in 1:length(formulas)) {
    
    parametros[i] <- strsplit(formulas[i], split = c("=~"))[[1]][2]
    
    parametros[i] <- strsplit(parametros[i][[1]], split = c("\\+"))
    for( j in 1:length(parametros[i][[1]])) {
      parametros[i][[1]][j] <- trimws(parametros[i][[1]][j])
    }
    ### 
    difference  <- unique(unlist(c(difference, setdiff(pegar5[[i]], parametros[i][[1]]))))
    
  }
  return(difference)
} 
retirarVarPadMenor  <- function(i) {
  v_l_aux<-NULL
  rhs_aux<-NULL
  v_l_aux<-retirar[i,1]
  rhs_aux<-retirar[i,2]
  
  k<-grep(rhs_aux, varia[[v_l_aux]])
  varia[[v_l_aux]][k]<-NA
  varia[[v_l_aux]]<-varia[[v_l_aux]][complete.cases(varia[[v_l_aux]])]
  return(varia)
}
tempos2 <- list()
indiceTempos2 <- 1
imprimirTempoPonderacao <- function(tempo_ponderar_inicio, msg) {
  tempo_ponderar_final<-NULL
  tempo_ponderar_final<-Sys.time()
  duracao_ponderacao<-NULL
  duracao_ponderacao<-tempo_ponderar_final-tempo_ponderar_inicio
  print(msg)
  print(duracao_ponderacao)
  tempos2$time[indiceTempos2] <- duracao_ponderacao
  tempos2$message[indiceTempos2] <- msg
  indiceTempos2 <- indiceTempos2 + 1
  return (tempos2)
}
varExcluidas <- function(varargs) {
  var_nao_queremos_inputar<-varargs
  variaveis_do_modelo<-gsub(var_nao_queremos_inputar
                            ,NA, variaveis_do_modelo)
  variaveis_do_modelo<-variaveis_do_modelo[complete.cases(variaveis_do_modelo)]
  return(variaveis_do_modelo)
}
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
    print("inicializando estimator")
    estimator <- match.arg(estimator) 
    if(estimator=="ML") warning("Estimator 'ML' will not correct standard errors and chi-square statistic.")
    print("inicializando teste estimator")
    estimator.gamma <- match.arg(estimator.gamma) # Smoothing Gamma or not
    print("inicializando teste lavaan names")
    # Names of the observed variables (same for each group)
    ov.names <- lavaanNames(lavaan.fit, type="ov", group=1)
    print("inicializando dplus (matrix duplication)")
    # The MP-inverse duplication matrix is handy for removing redundancy
    Dplus <- lavaan::lav_matrix_duplication_ginv( length(ov.names) )
    print("inicializando criacao formula")
    # Create a formula that includes all observed variables for svymean
    ov.formula <- as.formula(paste("~", paste(ov.names, collapse="+")))
    print("inicializando inspecao ngroups")
    # <no. group>-sized lists that will contain the asy. covariance matrix,
    #  and sample covariance matrix and mean vector
    ngroups <- lavInspect(lavaan.fit, "ngroups")
    print("inicializando Gamma vector")
    Gamma <- vector("list", ngroups)
    print("inicializando sample.cov")
    sample.cov <- vector("list", ngroups)
    print("inicializando sample.mean")
    sample.mean <- vector("list", ngroups)
    print("inicializando sample.nobs")
    sample.nobs <- vector("list", ngroups)
    
    print("inicializando for ngroups")
    for(g in seq(ngroups)) {
      if(ngroups > 1) {
        # Use survey::subset to create data groups
        print("inicializando survey.design.g")
        survey.design.g <- 
          subset(survey.design, eval(parse(text=sprintf("%s == '%s'", 
                                                        lavInspect(lavaan.fit, "group"), 
                                                        lavInspect(lavaan.fit, "group.label")[[g]]))))
      } 
      else { # In case of no groups, just use the original survey design object.
        print("inicializando survey.design.g no groups")
        survey.design.g <- survey.design  
      }
      
      # Function that takes survey design and returns the Gamma & observed moments
      get.stats.design <- function(survey.design.g, sample.nobs.g) {
        print("inicializando matrix transformation")
        ##essa operacao e demorada
        sample.cov.g <- as.matrix(svyvar(ov.formula, design=survey.design.g, na.rm=TRUE))  
        print("inicializando covariances matrix as attribute")
        # survey package returns the variance matrix of the (co)variances as attr:
        Gamma.cov.g <- attr(sample.cov.g, "var")
        print("inicializando remocao covariances")
        ##operacao pode ser demorada
        # Remove (co)variances wrt redundant elements of covs; not used by lavaan. 
        Gamma.cov.g <- Dplus %*% Gamma.cov.g %*% t(Dplus)
        print("inicializando remocao mean vector")
        # Same for mean vector
        sample.mean.g <- svymean(ov.formula, design=survey.design.g, na.rm=TRUE)  
        Gamma.mean.g <- attr(sample.mean.g, "var")
        # Join asy. variance matrices for means and covariances
        # TODO add offdiag
        print("inicializando lavaan matrix")
        Gamma.g <- lavaan::lav_matrix_bdiag(Gamma.mean.g, Gamma.cov.g)
        Gamma.g <- Gamma.g * sample.nobs.g # lavaan wants nobs * Gamma.
        # Since the above nonparametric estimate of Gamma can be unstable, Yuan
        # and Bentler suggested a model-smoothed estimate of it, optional here:
        if(estimator.gamma == "Yuan-Bentler") {
          print("inicializando yuan-bentler residuals")
          r <- get.residuals(lavaan.fit) # Iff these asy = 0, all will be well...
          Gamma.g <- Gamma.g + (sample.nobs.g/(sample.nobs.g - 1)) * (r %*% t(r))
        }
        # Get rid of attribute, preventing errors caused by lazy evaluation
        # (This has to be at the end or lazy evaluation mayhem will ensue)
        attr(sample.cov.g, "var") <- NULL
        tmp  <- as.vector(sample.mean.g)
        names(tmp) <- names(sample.mean.g)
        sample.mean.g <- tmp
        print("inicializando list gamma")
        list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
      }
      # The data may be a list of multiply imputed datasets
      if(!any(class(survey.design.g) == "svyimputationList")) {
        print("inicializando lavInspect lavaan.fit")
        # If no imputations, just use usual no. observations and asy variance
        sample.nobs.g <- lavInspect(lavaan.fit, "nobs")[[g]] 
        print("inicializando get.stats")
        stats <- get.stats.design(survey.design.g, sample.nobs.g)
        
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
        print("inicializando gamma reduction/calcs")
        # Variance estimates depend on within- and between-imputation variance:
        Gamma.within  <- Reduce(`+`, lapply(stats.list, `[[`, 'Gamma.g')) / m
        Gamma.between <- cov(cbind(mean.df, cov.df))
        Gamma.g <- Gamma.within + ((m + 1)/m) * Gamma.between
        
        # set stats with multiple imputation point and variance estimates
        stats <- list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
      }
      print("inicializando juncao lists")
      # Augment the list for this group
      Gamma[[g]] <- stats$Gamma.g
      sample.cov[[g]] <- stats$sample.cov.g
      sample.mean[[g]] <- stats$sample.mean.g
      sample.nobs[[g]] <- sample.nobs.g
      stats <- NULL
    } # End of loop over groups
    print("end of loop ngroup")
    Dplus <- NULL
    print("preenchendo new.call")
    new.call <- lavInspect(lavaan.fit, "call")
    new.call$data <- NULL                # Remove any data argument
    new.call$sample.cov <- sample.cov    # Set survey covariances
    new.call$sample.mean <- sample.mean  # Set survey means
    new.call$sample.nobs <- sample.nobs  
    new.call$estimator <- estimator  # Always use Satorra-Bentler or WLS estimator
    
    if(substr(estimator, 1, 2) == "ML") { # ML, robust options
      # Set asymptotic covariance matrix of sample means and covariances
      new.call$NACOV <- Gamma  
    }
    if(estimator %in% c("WLS", "DWLS")) {
      # Weighted Least Squares, adjust the weight matrix: MP inverse of Gamma
      # Note that Gamma may be singular.
      new.call$WLS.V <- lapply(Gamma, ginv)
    }
    print("preenchendo new.fit")
    ##outra operacao demorada
    new.fit <- eval(as.call(new.call)) # Run lavaan with the new arguments
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
    print("finalizando new.fit")
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
########################################################################################
###################################END of function######################################
########################################################################################
TempoTotal<-Sys.time()
########################################################################################
#### ATUALIZAR O PREENCHIMENTO    ######################################################

Inicio_Sintaxe_SEM<-Sys.time()

###inicialização variáveis para doParallel
no_cores <- (detectCores() -1)


# baixar banco de dados sem labels

data_raw <- read.spss("FPB_Student%20Satisfaction%202017_DB_V1_Gradua%C3%A7%C3%A3o.sav",
                      to.data.frame=T,
                      header=T,
                      use.value.labels = F,
                      use.missings = T)

## baixar os dados com labels:

data_value_labels <- read.spss("FPB_Student%20Satisfaction%202017_DB_V1_Gradua%C3%A7%C3%A3o.sav",
                               to.data.frame=T,
                               header=T,
                               use.value.labels = T,
                               use.missings = T)


##   diretorio para onde exportar os resultados

export<-"/home/gabrielpehls/sem/outputs"

##   quais quebras usar

quebras<-c("VERTICAL", "TYPE", "AnoIngresso","Finished")


### modelo descrito completo

model_completo <- '
campus=~INF_01_INF_01_01+INF_01_INF_01_02+INF_01_INF_01_04+INF_01_INF_01_06+INF_01_INF_01_08+INF_01_INF_01_09+INF_01_INF_01_11+INF_01_INF_01_13+INF_01_INF_01_15+INF_01_INF_01_12+INF_01_INF_01_19+INF_01_INF_01_21
classroom=~INF_02_INF_02_03+INF_02_INF_02_04+INF_02_INF_02_05+INF_02_INF_02_06+INF_02_INF_02_07
it_labs=~INF_03_INF_03_01+INF_03_INF_03_02+INF_03_INF_03_03+INF_03_INF_03_04+INF_03_INF_03_06+INF_03_INF_03_07
specific_labs=~INF_04_INF_04_01+INF_04_INF_04_02+INF_04_INF_04_03+INF_04_INF_04_06
library=~INF_05_INF_05_02+INF_05_INF_05_04+INF_05_INF_05_10+INF_05_INF_05_09+INF_05_INF_05_07
blackboard=~INF_06_INF_06_01+INF_06_INF_06_02+INF_06_INF_06_04+INF_06_INF_06_07+INF_06_INF_06_08
faculty=~PROF_02_PROF_02_01+PROF_02_PROF_02_02+PROF_02_PROF_02_03+PROF_02_PROF_02_04+PROF_02_PROF_02_05+PROF_02_PROF_02_06+PROF_02_PROF_02_07+PROF_02_PROF_02_08+PROF_02_PROF_02_20
coord=~COORD_02_COORD_02_01+COORD_02_COORD_02_02+COORD_02_COORD_02_03+COORD_02_COORD_02_04+COORD_02_COORD_02_20
call=~CALL_01_CALL_01_01+CALL_01_CALL_01_02+CALL_01_CALL_01_03+CALL_01_CALL_01_04
employ=~EMPL_02_EMPL_02_01+EMPL_02_EMPL_02_02+EMPL_02_EMPL_02_03+EMPL_02_EMPL_02_04+EMPL_02_EMPL_02_05+EMPL_02_EMPL_02_8
financial_serv=~FSERV_01_FSERV_01_01+FSERV_01_FSERV_01_02+FSERV_01_FSERV_01_03+FSERV_01_FSERV_01_04
program=~PROG_02_PROG_02_01+PROG_02_PROG_02_02+PROG_02_PROG_02_03+PROG_02_PROG_02_04+PROG_02_PROG_02_05
disc_online=~ELEAR_02_ELEAR_02_01+ELEAR_02_ELEAR_02_02+ELEAR_02_ELEAR_02_03+ELEAR_02_ELEAR_02_04+ELEAR_02_ELEAR_02_05
comunic=~COMM_01_COMM_01_01+COMM_01_COMM_01_01.0+COMM_01_COMM_01_01.1+COMM_01_COMM_01_01.2+COMM_01_COMM_01_02+COMM_01_COMM_01_03+COMM_01_COMM_01_06+COMM_01_COMM_01_08
std_services=~SSERV_01_SSERV_01_01+SSERV_01_SSERV_01_02+SSERV_01_SSERV_01_03+SSERV_01_SSERV_01_04
inter=~INTL_01_INTL_01_01+INTL_01_INTL_01_02+INTL_01_INTL_01_03+INTL_01_INTL_01_04+INTL_01_INTL_01_05+INTL_01_INTL_01_06+INTL_01_INTL_01_07
image=~IMG_01_IMG_01_01+IMG_01_IMG_01_02+IMG_01_IMG_01_03
#
infra=~campus+classroom+it_labs+specific_labs+library+blackboard
SAT_00_02~infra+image+inter+std_services+comunic+program+financial_serv+employ+call+coord+disc_online
SAT_00_01~infra+image+inter+std_services+comunic+program+financial_serv+employ+call+coord+disc_online+SAT_00_02
NPS_01_01~SAT_00_01+image'

constructs <- c("campus","classroom","it_labs","specific_labs","library","blackboard","disc_online")
var_dep_eq<-c("SAT_00_02","SAT_00_01","NPS_01_01")
segunda_ordem<-c("infra")


soh_val<-NULL
soh_val<-sapply(strsplit(model_completo, split= "\\#"),
                function(x) x[length(x)])
print(lev)
print(model_completo)

######################## FIM DA ATUALIZACAO DO PREENCHIMENTO #####################################
##################################################################################################


# Ponderar os resultados
svy.df<-svydesign(id=~ResponseId, 
                  weights=~Weight_BRZ,
                  data=data_raw)

#pegar as variable labels
var_labels<-data.frame(attr(data_raw, 'variable.labels'))
var_labels<-cbind(var_labels, rownames(var_labels))
rownames(var_labels)<-NULL
colnames(var_labels)<-c("label","codigo")

#################################################################
## Vamos imputar valores para os missing values      ############


pegar1<-NULL
pegar1<-unlist(strsplit(model_completo,
                        split=c("\n")))
pegar2<-NULL
pegar2<-grep("=~", pegar1, value=T)
#### vamos tirar a ultima formula ("infra=~campus+classroom+it_labs+specific_labs+library+blackboard" ) 
#### que não pertence ao conjunto de constructos
pegar2 <- pegar2[-length(pegar2)]


# agora vamos dividir entre constuctos e 'preditores'

pegar3<-list()
pegar4<-c()


for(i in 1:length(pegar2)){
  pegar3[i]<-strsplit(pegar2[i],
                      split=c("=~"))
  pegar4[i]<-trimws(pegar2[[i]][1])
  pegar4<-trimws(pegar4)
}

# separar os preditores dentro de cada constructo

pegar5<-list()

for(i in 1:length(pegar3)){
  
  ## vamos dividir em cada elemento da parte dos preditores 
  
  pegar5[i]<-strsplit(pegar3[[i]][2],
                      split=c("\\+"))
  
}

lista_variaveis<-NULL
lista_variaveis<-unlist(pegar5)

### comeca a verificar se as varivaveis observadas do modelo estao no banco de dados


variaveis_lat<- setdiff(lista_variaveis, colnames(data_value_labels))
variaveis_lat<- setdiff(variaveis_lat, constructs)

if (length(variaveis_lat) > 0) {
  latent<-NULL
  latent<-unlist(strsplit(model_completo,
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
  
  print("dividir entre constuctos e preditores:ok")
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
  ### NO TESTE DE MATRIZ DE COVARIANCIAS
  k<-NULL
  var_aux<-NULL
  res<-NULL
  aux<-NULL
  
  for(i in 1:length(variaveis_lat)){
    
    variaveis_lat_aux<-variaveis_lat[i]
    res <- lapply(varia, function(ch) grep(variaveis_lat_aux, ch))
    aux<-sapply(res, function(x) length(x) > 0)
    
    if(sum(aux)>0){
      
      var_aux<-names(which(aux==T))
      k<-grep(variaveis_lat_aux, varia[[var_aux]])
      varia[[var_aux]][k]<-NA
      varia[[var_aux]]<-varia[[var_aux]][complete.cases(varia[[var_aux]])]
      
    }
    # variaveis_lat_aux<-NULL
    k<-NULL
    var_aux<-NULL
    res<-NULL
    aux<-NULL
  }
  
  model_completo <- NULL
  for(i in 1:length(varia)){
    
    va_lt<-NULL
    
    va_lt[i]<-paste(names(varia)[i],"=~")
    
    modelo<-NULL
    var_lat<-(paste(unlist(varia[[i]]),"+"))
    modelo<-paste(var_lat, collapse = " ")
    modelo<-substr(modelo, 1, nchar(modelo)-2)
    modelo<-paste(va_lt[i], modelo,"\n")
    model_completo<-paste(model_completo, modelo)
    print(lev)
    print(model_completo)
  }
  ## dividir variáveis
  pegar1<-NULL
  pegar1<-unlist(strsplit(model_completo,
                          split=c("\n")))
  pegar2<-NULL
  pegar2<-grep("=~", pegar1, value=T)
  
  # agora vamos dividir entre constuctos e 'preditores'
  
  pegar3<-list()
  pegar4<-c()
  
  
  for(i in 1:length(pegar2)){
    pegar3[i]<-strsplit(pegar2[i],
                        split=c("=~"))
    pegar4[i]<-trimws(pegar2[[i]][1])
    pegar4<-trimws(pegar4)
  }
  
  # separar os preditores dentro de cada constructo
  
  pegar5<-list()
  
  for(i in 1:length(pegar3)){
    
    ## vamos dividir em cada elemento da parte dos preditores 
    
    pegar5[i]<-strsplit(pegar3[[i]][2],
                        split=c("\\+"))
    
  }
  
  lista_variaveis<-NULL
  lista_variaveis<-unlist(pegar5)
  lista_variaveis <- trimws(lista_variaveis)
  model_completo <- paste(model_completo, soh_val,sep = " #")
}

variaveis_do_modelo<-Reduce(intersect, list(lista_variaveis,colnames(data_raw)))


## transformar as vaiaveis usadas em numericas
indices_colunas<-NULL
indices_colunas<-which(colnames(data_raw) %in% variaveis_do_modelo)

for(j in 1:length(indices_colunas)){
  j<-indices_colunas[i]
  data_raw[,j] <-  as.numeric(as.character(data_raw[,j]))
}
print(lev)
print(model_completo)
#######################################################################

### ATENCAO: CASO NAO QUEIRAMOS INPUTAR DADOS PARA ALGUMA VARIAVEL ###

######################################################################

## USAR ABAIXO:
## QUAIS

#variaveis_do_modelo<-varExcluidas(c('',''))

############################################################################


########################################################
#### ANTES: CALCULAR O T2B PARA CADA VARIAVEL   ########
#######################################################

indices_colunas<-NULL
# indices_colunas<-which(colnames(data_raw) %in% variaveis_do_modelo)
indices_colunas<-variaveis_do_modelo


## Finalmente, vamos adicionar a media nos NA's das variaveis
## Colocaremos as medias em funcao das verticais


# imputar medias gerais pra cada variavel

data_raw$concatenar<-do.call(paste, data_raw[,quebras])

# percorrer cada individuo, em cada variavel, ver se esta nulo e colocar
# a media da variavel conforme a vertical que pertence

for(k in 1:length(indices_colunas)){
  
  j<-indices_colunas[k]
  
  teste_nan<-NULL
  
  teste_nan<-aggregate(data_raw[,j],
                       list(data_raw$concatenar), 
                       mean,
                       na.rm=T)
  
  
  
  svy.df<-NULL
  svy.df<-svydesign(id=~ResponseId, 
                    weights=~Weight_BRZ,
                    data=data_raw)
  
  
  media_geal_aux<-NULL
  media_geal_aux<-svymean(~data_raw[,j], 
                          design=svy.df,
                          na.rm=T,
                          deff=T)[1]
  
  teste_nan$x[is.nan(teste_nan$x)]<-media_geal_aux
  
  # 
  # 
  # if(sum(is.nan(teste_nan$x))>0){
  #   
  #           concat_NoN_NaN<-NULL
  #           concat_NoN_NaN<-teste_nan[complete.cases(teste_nan$x),1]
  #           data_raw_aux<-NULL
  #           data_raw_aux<-data_raw[data_raw$concatenar %in% concat_NoN_NaN,]
  #           
  #           svy.df<-svydesign(id=~ResponseId, 
  #                             weights=~Weight_BRZ,
  #                             data=data_raw_aux)
  #           medias_aux<-NULL
  #           medias_aux<-svyby(~eval(as.symbol(j)),
  #                             ~concatenar,
  #                             svy.df,
  #                             svymean, 
  #                             na.rm =T)  
  #           
  # }else{
  #   
  #   svy.df<-svydesign(id=~ResponseId, 
  #                     weights=~Weight_BRZ,
  #                     data=data_raw)
  #   medias_aux<-NULL
  #   medias_aux<-svyby(~eval(as.symbol(j)),
  #                     ~concatenar,
  #                     svy.df,
  #                     svymean, 
  #                     na.rm =T)  
  # }
  
  
  for(i in 1:nrow(data_raw)){
    
    vertical_aux<-data_raw$concatenar[i]
    
    if(is.na(data_raw[i,j])){
      # data_aux<-NULL
      # data_aux<-data_raw[data_raw$concatenar==vertical_aux,]
      # 
      ### se toda a variavel estiver vazia na vertical,
      ### entao temos que imputar a media geral
      
      data_raw[i,j]<-teste_nan$x[teste_nan$Group.1==vertical_aux]
      
      # if(sum(is.na(data_raw[data_raw$concatenar==vertical_aux,j]))==length(data_aux[,j])){
      #   
      #   svy.df<-NULL
      #   svy.df<-svydesign(id=~ResponseId, 
      #                     weights=~Weight_BRZ,
      #                     data=data_raw)
      #   
      #   data_raw[i,j]<- svymean(data_raw[,j], design=svy.df, na.rm=T)[1]
      #   
    }
  }
}



# Atualizar a ponderacao os resultados
svy.df<-svydesign(id=~ResponseId, 
                  weights=~Weight_BRZ,
                  data=data_raw)

###################################################################################################
###                                                                     ###########################
#### A partir daqui fica automatico                                     ###########################
##                                                                      ###########################
##                                                #################################################


### A gente vai fazer para cada quebra necessaria    #############################################

for(q in 1:length(quebras)){
  ## quebra 'q'
  #q<-1
  
  quebra<-quebras[q]
  niv<-levels(as.factor(data_value_labels[,quebras[q]]))
  
  
  for(lev in 1:length(niv)){  
    #lev <- 1
    ###  nivel 'niv' da quebra 'q'
    nivel<-levels(as.factor(data_raw[,quebras[q]]))[lev]
    ## Nome do nivel da quebra
    nivel_label<-levels(as.factor(data_value_labels[,quebras[q]]))[lev]
    
    if (is.na(nivel)) break
    
    data<-NULL
    data<-data_raw[data_raw[,quebra]==nivel,]
    data<-data[,c(var_dep_eq, "ResponseId",
                  variaveis_do_modelo,"Weight_BRZ")]
    
    
    n_obs<-nrow(data)
    if (n_obs < 50) {
      print(paste(paste(nivel_label, " Nobs="),n_obs))
      break
    }
    print(nivel_label)
    
    svy.df<-NULL
    svy.df<-svydesign(id=~ResponseId, 
                      weights=~Weight_BRZ,
                      data=data)
    
    
    ############################################################################
    ######### BANCO AUXILIAR PARA O TOP 2 BOX      #############################
    
    
    # ### banco auxiliar para top 2 box
    # 
    data_aux_top<-data_raw[data_raw[,quebra]==nivel,] %>%
      select(indices_colunas, "Weight_BRZ", "ResponseId", quebras)
    # executado na linha 136
    
    
    for(i in 1:length(indices_colunas)){
      data_aux_top[,i]<- ifelse(data_aux_top[,i]>3,"T2B", 
                                ifelse(data_aux_top[,i]<3,"B2B","Neutral"))
      
    }
    
    
    i<-NULL
    # Ponderar os resultados
    svy.df.auxtop<-NULL
    svy.df.auxtop<-svydesign(id=~ResponseId, 
                             weights=~Weight_BRZ,
                             data=data_aux_top)
    df.top2box<-data.frame()
    
    for(i in 1:length(indices_colunas)){
      # i<-52
      nome_aux_t2b<-NULL    
      nome_aux_t2b<-colnames(data_aux_top[i])
      aux_t2box<-NULL
      aux_t2box<-as.data.frame(prop.table(svytable(~eval(as.symbol(colnames(data_aux_top)[i])),
                                                   design = svy.df.auxtop))*100)
      
      df.top2box_aux<-NULL
      df.top2box_aux<-cbind(nome_aux_t2b,aux_t2box$Freq[aux_t2box[,1]=="T2B"])
      df.top2box<-rbind(df.top2box,df.top2box_aux)
    }
    
    colnames(df.top2box)<-c("variaveis", "T2B_result")
    
    ###############################################################################
    ###############################################################################
    
    #### o simbolo "=~" significa que a a variavel latente y se manifesta pelas variaveis
    #### observadas x1, x2,...xk
    
    
    HS.model<-NULL
    HS.model <- substr(model_completo,1,
                       regexpr("#",model_completo)-1)
    
    
    #########################################################################
    ### VAMOS TESTAR SE A MATRIZ DE COVARIANCIA EH DEFINIDA POSITIVA #######
    
    
    teste_covariancia<-(is.positive.definite(cov(data[,indices_colunas])) 
                        & n_obs>length(indices_colunas))
    
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
      
      correl_tirar_var<-cor(data[,indices_colunas])
      
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
      
      
      if(t2b_var2$T2B_result.x==t2b_var2$T2B_result.y){
        
        sorteio_aux<-NULL
        set.seed(1234)
        sorteio_aux<-sample(c(1,2),1)
        t2b_var2$retirar<-as.character(t2b_var2[1,sorteio_aux])
        
      }else{
        
        
        t2b_var2$retirar<-ifelse(as.numeric(as.character(t2b_var2$T2B_result.x))<as.numeric(as.character(t2b_var2$T2B_result.y)), 
                                 t2b_var2$retirar<-as.character(t2b_var2$Var2),
                                 t2b_var2$retirar<-as.character(t2b_var2$Var1))
        
        
        
      }
      ## coloca num vetor quais variaveis ficaram de fora
      retirar_var_aux<-NULL
      retirar_var_aux<-t2b_var2$retirar
      variaveis_PROBLEMA_MATRIZ[ind]<-retirar_var_aux
      
      coord_retirar<-NULL
      coord_retirar<-which(indices_colunas %in% retirar_var_aux)
      indices_colunas<-indices_colunas[-coord_retirar]
      
      
      
      ## testa se eh positiva definida
      teste_covariancia<-(is.positive.definite(cov(data[,indices_colunas])) 
                          & n_obs>length(indices_colunas))
      
      if(teste_covariancia==F | length(variaveis_PROBLEMA_MATRIZ)>0){
        ind<-ind+1
      }else{
        ind<-ind
      }
    }
    
    
    
    ###############################################
    
    ## TIRAMOS AS VARIAVEIS QUE COM A MULTICOLINEARIDADE
    
    ## AFETAVAM NA MATRIZ.
    
    ##  AGORA PRECISAMOS TIRAR AS VARIAVEIS QUE SAIRAM DO MODELO INICIAL
    ### ISSO SO VAI ACONTECER SE FORAM FEITAS ALTERACOES EM FUNCAO
    ## DO TAMANHO DE AMOSTRA OU DE A MATRIZ NAO SER DEFINIDA POSITIVA
    
    if(ind>1){
      
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
      
      print("dividir entre constuctos e preditores:ok")
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
      ### NO TESTE DE MATRIZ DE COVARIANCIAS
      variaveis_PROBLEMA_MATRIZ_aux<-NULL
      k<-NULL
      var_aux<-NULL
      res<-NULL
      aux<-NULL
      
      for(i in 1:length(variaveis_PROBLEMA_MATRIZ)){
        
        variaveis_PROBLEMA_MATRIZ_aux<-variaveis_PROBLEMA_MATRIZ[i]
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
    
    fit_aux<-NULL
    
    fit_aux <- cfa(HS.model,
                   data=data,
                   likelihood="wishart", std.lv=T)
    
    
    fit <- NULL
    tempo_ponderar_inicio<-Sys.time()
    
    modelo_partial<-list()
    modelo_name<-list()
    ###########################################################################
    print("inicializando lavaan survey")
    modelo_partial[[1]]<- lavaan.survey(fit_aux, 
                                        svy.df) 
    print("finalizando lavaan survey")
    # 
    # fit <- lavaan.survey(fit_aux, 
    #                    svy.df) 
    ######################################
    
    tempos2 <- imprimirTempoPonderacao(tempo_ponderar_inicio, "Modelo parcial: ")
    indiceTempos2 <- indiceTempos2 +1
    ######################################
    # tempo_ponderar_final<-Sys.time()
    # 
    # 
    # duracao_ponderacao<-tempo_ponderar_final-tempo_ponderar_inicio
    # print("A ponderação do modelo parcial é "+duracao_ponderacao)
    # 
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
    # summ_cfa_modelo_comp<-summary(fit, 
    #                               fit.measures=T,
    #                               rsquare=T, standardized=T)
    # 
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
    
    
    
    if( fit_modelo_partial["rmsea.ci.upper"]<0.08 & fit_modelo_partial["srmr"]<0.08 & fit_modelo_partial["cfi"]>0.8 & nrow(est_aux)<1 & r2_aux<1 & var_aux<1 & est_pad_aux<1 & chisq_df<5){
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
    #modelo_partial<-list()
    modelo_name<-list()
    
    # modelo_partial[[1]]<-fit
    modelo_name[[z]]<-HS.model
    
    print("vai entrar no while da CFA")
    
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
      
      print("dividir entre constuctos e preditores:ok")
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
      
      print("lista varia")
      
      # a lista 'varia' tem o nome dos constructos entao
      # vamos nos dedicar a tirar ou colocar variaveis 
      # com objetivo de melhorar o ajuste
      
      # tirar as variaveis nao significativas
      print("")
      estimativas<-parameterEstimates(modelo_partial[[z]])[complete.cases(parameterEstimates(modelo_partial[[z]])),]
      est_aux<-NULL
      est_aux<-estimativas[(estimativas$pvalue>0.05
                            & estimativas$op=="=~")| (estimativas$est<0 &
                                                        estimativas$op=="=~")  ,]
      est_aux<-est_aux[complete.cases(est_aux$z),]
      
      retirar<-NULL
      retirar<-cbind(est_aux$lhs, est_aux$rhs)
      i <- 1
      ###############################
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
      
      ###############################
      
      print("retirar variaveis nao significativas:ok")
      estimativas2<-NULL
      #### Retirar variaveis cuja padronizacao for menor que 0.25
      estimativas2<-parameterEstimates(modelo_partial[[z]], standardized = T)[complete.cases(parameterEstimates(modelo_partial[[z]])),]
      
      est_aux2<-estimativas2[(estimativas2$std.lv<0.25
                              & estimativas$op=="=~"),]
      retirar<-NULL
      retirar<-cbind(est_aux2$lhs, est_aux2$rhs)
      i <- 1
      ###############################
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
      
      ###############################
      
      print("retirar variaveis cujo pad for menor que 0.25:ok")
      
      ## alteracoes para incluir variaveis em constructos
      MI<-NULL
      MI<-modificationIndices(modelo_partial[[z]], sort. = T)
      # adc_aux<-subset(MI, mi>100 & op=="=~" & epc>0)
      # adc_aux <- adc_aux[order(adc_aux[,'rhs'],-adc_aux[,'epc']),]
      # adc_aux <- adc_aux[!duplicated(adc_aux$rhs),]
      # adicionar<-cbind(adc_aux$lhs, adc_aux$rhs)
      # 
      # ## Esses indices podem ser vistos como a estatistica chisq com um 
      # ## grau de liberdade. Para cada parametro especficiado, fixo,
      # ## um valor eh dado, e representa a queda esperada no valor do chisq
      # ## geral se o parametro nao fosse mais fxo e sim livremente estimavel,
      # ## numa proxima vez em que o modelo fosse executado
      # 
      # ## primeiro, retirar onde estao por enquanto
      # if(nrow(adicionar)>0){
      #   for(i in 1:nrow(adicionar)){
      #     
      #     v_l_aux<-NULL
      #     rhs_aux<-NULL
      #     k<-NULL
      #     rhs_aux<-adicionar[i,2]
      #     
      #     tab_aux<-lapply(varia, function(ch) grep(rhs_aux, ch))
      #     
      #     # verificar se a variavel nao esta em nenhum constructo
      #     contar<-0
      #     soma<-NULL
      #     
      #     for(k in 1:length(tab_aux)){
      #       soma<-length(tab_aux[[k]])
      #       contar<-contar+soma
      #     }
      #     
      #     if(contar>0){
      #       
      #       v_l_aux<-names(which(tab_aux>0))
      #       
      #       k<-grep(rhs_aux, varia[[v_l_aux]])
      #       varia[[v_l_aux]][k]<-NA
      #       varia[[v_l_aux]]<-varia[[v_l_aux]][complete.cases(varia[[v_l_aux]])]
      #     }
      #   }
      #   print("retirou as variaveis de onde estavam:ok")
      #   ## colocar em outras variaveis latentes
      #   
      #   for(i in 1:nrow(adicionar)){
      #     
      #     v_l_aux<-NULL
      #     rhs_aux<-NULL
      #     k<-NULL
      #     v_l_aux<-adicionar[i,1]
      #     rhs_aux<-adicionar[i,2]
      #     
      #     k<-length(varia[[v_l_aux]])+1
      #     varia[[v_l_aux]][k]<-rhs_aux
      #   }
      # }
      # 
      # print("Colocar variaveis em outros constructos:ok")
      # 
      
      ### Ainda o Modification Indices pode nos ajudar a identificar
      ### problemas de multicolinearidade. 
      
      ## Pegar os tres casos que evidenciaram melhorias ao adicionar covariancias 
      multi_aux<-head(subset(MI, mi>130 & op=="~~"),3)
      
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
        # estimativas_aux2.est.x se refere multi_aux.lhs
        # estimativas_aux2.est.y se refere multi_aux.rhs
        
        
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
      ####
      
      print("Tirar multicolinearidade:ok")
      
      ### Temos que ver ainda se alguma variavel latente nao ficou vazia:
      
      varia<-Filter(length, varia)
      
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
      
      cfa_GSP_aux<-NULL
      cfa_GSP_aux<-cfa(model = modelofinal,
                       likelihood="wishart",std.lv=T,
                       data=data)
      
      tempo_ponderar_inicio<-NULL
      tempo_ponderar_inicio<-Sys.time()
      
      z<-z+1  
      modelo_partial[[z]]<-lavaan.survey(cfa_GSP_aux, 
                                         svy.df) 
      
      # cfa_GSP<-NULL
      # cfa_GSP <- lavaan.survey(cfa_GSP_aux, 
      #                          svy.df) 
      ######################################
      
      tempos2 <- imprimirTempoPonderacao(tempo_ponderar_inicio, paste(paste("Modelo parcial ",z),": "))
      indiceTempos2 <- indiceTempos2 +1
      ######################################
      # tempo_ponderar_final<-NULL
      # tempo_ponderar_final<-Sys.time()
      # duracao_ponderacao<-NULL
      # duracao_ponderacao<-tempo_ponderar_final-tempo_ponderar_inicio
      # print("A duracao da ponderacao do modelo parcial "+z+" e "+duracao_ponderacao)
      # 
      bic[z]<-fitMeasures(modelo_partial[[z]], "bic")
      print(bic[z])
      
      #z<-z+1  
      #modelo_partial[[z]]<-cfa_GSP
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
      ### verificar se tem alguma estimativa padronizada menor que 0.25
      estimativas_pad<-NULL
      est_pad_aux<-NULL
      chisq_df<-NULL
      estimativas_pad<-parameterEstimates(modelo_partial[[z]],standardized = T)[complete.cases(parameterEstimates(modelo_partial[[z]])),]
      est_pad_aux<-sum(estimativas_pad$std.all[estimativas_pad$op=="=~"]<0.25, na.rm=T)
      fit_modelo_partial <- NULL
      fit_modelo_partial <- fitMeasures(modelo_partial[[z]])
      chisq_df<-fit_modelo_partial["chisq"]/fit_modelo_partial["df"]
      fit<-NULL
      # fit<-cfa_GSP
      
      parar<-0
      if(z>2){
        parar<-sum(bic[z-1]==bic[z-2])
      }else{
        parar<-0
      }
      
      if( (fit_modelo_partial["rmsea.ci.upper"]<0.08 & nrow(est_aux)<1 & r2_aux<1 & var_aux<1 & est_pad_aux<1 & chisq_df<5) |  (z>30 | parar==1) ){
        rodar<-0
      }else{
        rodar<-1
      }
      adc_aux<-NULL
      aval_multi1<-NULL
      aval_multi2<-NULL
      
      print("Ainda nao saiu do while")
      print(lev)
      print(soh_val)
      
    }
    
    
    print("saiu do while da CFA")
    
    if(z>1){
      for(i in 2:(z)){
        fit_modelo_partial<- NULL
        fit_modelo_partial<-fitMeasures(modelo_partial[[i]])
        alfa_cronbach[i]<-reliability(modelo_partial[[i]])[1,ncol(reliability(modelo_partial[[i]]))]
        razao_chisq[i]<-fit_modelo_partial["chisq"]/fit_modelo_partial["df"]
        aic[i]<-fit_modelo_partial["aic"]
        rmsea[i]<-fit_modelo_partial["rmsea"]
        cfi[i]<-fit_modelo_partial["cfi"]
        ave_min[i]<-min(reliability(modelo_partial[[i]])[5,])
        ## matriz de validade discriminante
        mat_validDis<-inspect(modelo_partial[[i]], "cor.lv")
        diag(mat_validDis)<-NA
        valid_discr[i]<-sum(mat_validDis>0.95, na.rm=T)
      }
      
      
      ###pegar o modelo com maior cfi
      
      melhor<-which.max(cfi)
      ## precisamos tirar os constructos que ficaram de fora:
      ## constructos que sobraram:
    }else{
      melhor<-1
    }
    ## funcao que identifica diferenca entre dois vetores
    outersect <- function(x, y) {
      sort(c(setdiff(x, y),
             setdiff(y, x)))
    }
    # quais constructos ficaram de fora
    
    diferentes<-outersect(nome_constructo,valat_originais)
    
    ## soh para o caso de detectar mais de 0 constructos que ficaram de fora
    
    print(" saiu do while e vai ver se tem constructos diferentes")
    print(lev)
    print(soh_val)
    
    if(length(diferentes)>0){
      ## tirando o constructo
      tirar_val<-strsplit(soh_val, split=diferente, fixed=F)
      
      verificar<-substring(tirar_val[[1]][2], 1, 1)
      ### se o primeiro character da segunda parte for '+',
      ## entao precisamos retirar
      if(verificar=="+"){
        tirar_val[[1]][2]<-sub('.', '', tirar_val[[1]][2])
      }
      soh_val<-paste(tirar_val[[1]],collapse="")
      
    }
    
    
    modelo_completo2<-paste(modelo_name[[melhor]], soh_val, sep = "#")
    
    print(" saiu do while e vai ver se tem constructos diferentes")
    print(lev)
    print(soh_val)
    
    
    rodar_modelo2_aux<-NULL
    rodar_modelo2_aux<-sem(modelo_completo2,  
                           likelihood="wishart",std.lv=T,
                           data=data)
    
    tempo_ponderar_inicio<-NULL
    tempo_ponderar_inicio<-Sys.time()
    
    survey.fit<-NULL
    survey.fit <- lavaan.survey(rodar_modelo2_aux, 
                                svy.df) 
    ######################################
    
    tempos2 <- imprimirTempoPonderacao(tempo_ponderar_inicio, "Survey.fit: ")
    indiceTempos2 <- indiceTempos2 +1
    ######################################
    # tempo_ponderar_final<-NULL
    # tempo_ponderar_final<-Sys.time()
    # duracao_ponderacao<-NULL
    # duracao_ponderacao<-tempo_ponderar_final-tempo_ponderar_inicio
    # print("A duracao da ponderacao da survey.fit e "+duracao_ponderacao)
    # 
    # summary(survey.fit, standardized=T,
    #         rsquare=T)
    # 
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
    constru_nao_sig<-regre_pad_aux[regre_pad_aux$op=="~" & (regre_pad_aux$pvalue>0.05 | regre_pad_aux$est<0),]
    constru_nao_sig<-cbind(constru_nao_sig$lhs, constru_nao_sig$rhs)
    
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
      
      print("primeira etapa da separacao variaveis")
      
      # separar os preditores dentro de cada constructo
      
      pegar_reg5<-list()
      var_dep_reg<-c()
      
      for(i in 1:length(pegar_reg3)){
        
        ## vamos dividir em cada elemento da parte dos preditores 
        
        pegar_reg5[i]<-strsplit(pegar_reg3[[i]][2],
                                split=c("\\+"))
        
        var_dep_reg[i]<-pegar_reg3[[i]][1]  
      }
      
      print("separar os preditores dentro de cada construco")
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
        
        va_lt<-NULL
        
        va_lt[i]<-paste(names(pegar_reg5)[i],"~")
        
        modelo<-NULL
        var_lat<-(paste(unlist(pegar_reg5[[i]]),"+"))
        modelo<-paste(var_lat, collapse = " ")
        modelo<-substr(modelo, 1, nchar(modelo)-2)
        modelo<-paste(va_lt[i], modelo,"\n")
        modelo_reg<-paste(modelo_reg, modelo)
      }
      
      
      
      modelo_completo3<-NULL
      modelo_completo3<-paste(modelo_name[[melhor]], modelo_reg)
      
      rodar_modelo3_aux<-sem(modelo_completo3,  
                             likelihood="wishart",std.lv=T,
                             data=data)
      
      tempo_ponderar_inicio<-NULL
      tempo_ponderar_inicio<-Sys.time()
      
      survey.fit2 <- lavaan.survey(rodar_modelo3_aux, 
                                   svy.df)
      ######################################
      
      tempos2 <- imprimirTempoPonderacao(tempo_ponderar_inicio, "Survey.fit2: ")
      indiceTempos2 <- indiceTempos2 +1
      ######################################
      # tempo_ponderar_final<-NULL
      # tempo_ponderar_final<-Sys.time()
      # duracao_ponderacao<-NULL
      # duracao_ponderacao<-tempo_ponderar_final-tempo_ponderar_inicio
      # print("A duracao da ponderacao da survey,fit2 e "+duracao)
      
      ## vamos pegar as estimacoes das regressoes
      regre_pad_aux<-NULL
      regre_pad_aux<-parameterEstimates(survey.fit2,
                                        standardized = T)[complete.cases(parameterEstimates(survey.fit)),]
      
      ## retornar as var depend e os respectivos constructos nao significatios
      constru_nao_sig<-NULL
      constru_nao_sig<-regre_pad_aux[regre_pad_aux$op=="~" & (regre_pad_aux$pvalue>0.05 | regre_pad_aux$est<0),]
      constru_nao_sig<-constru_nao_sig[complete.cases(constru_nao_sig),]  
      constru_nao_sig<-(cbind(constru_nao_sig$lhs, constru_nao_sig$rhs))
      
      
      
      ver_ajuste<-(nrow(constru_nao_sig)>0 & contador<5)
      qtos_prob[contador]<-nrow(constru_nao_sig)
      print(qtos_prob[contador])
      contador<-contador+1
      print(contador)
      
      
    }
    
    
    
    ### Modelo final
    tempo_ponderar_inicio<-NULL
    tempo_ponderar_inicio<-Sys.time()
    
    if(contador>1){
      rodar_modelo3_aux<-sem(modelo_completo3,  
                             likelihood="wishart",std.lv=T,
                             data=data)
      
      modelo_final_final<- lavaan.survey(rodar_modelo3_aux, 
                                         svy.df) 
    }else{
      rodar_modelo2_aux<-sem(modelo_completo2,  
                             likelihood="wishart",std.lv=T,
                             data=data)
      
      modelo_final_final<- lavaan.survey(rodar_modelo2_aux, 
                                         svy.df) 
      
      
    }
    ######################################
    
    tempos2 <- imprimirTempoPonderacao(tempo_ponderar_inicio, "Modelo final: ")
    indiceTempos2 <- indiceTempos2 +1
    
    ######################################
    # tempo_ponderar_final<-NULL
    # tempo_ponderar_final<-Sys.time()
    # duracao_ponderacao<-NULL
    # duracao_ponderacao<-tempo_ponderar_final-tempo_ponderar_inicio
    # print("A duracao da ponderacao final é "+ duracao_ponderacao)
    # summary(modelo_final_final,standardized=T,
    #         rsquare=T)
    # 
    
    
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
    est_modelo_final_val_aux$interp_pred<-NA
    est_modelo_final_val_aux$interp_pred[is.na(est_modelo_final_val_aux$interp_pred) & est_modelo_final_val_aux$op=="=~"]<-paste("O aumento de uma unidade na satisfacao da variavel latente",
                                                                                                                                 est_modelo_final_val_aux$label_lhs,"implica no aumento de",
                                                                                                                                 round(est_modelo_final_val_aux$est,4), "unidades na variavel observada",
                                                                                                                                 est_modelo_final_val_aux$label_rhs)
    
    est_modelo_final_reg_aux<-est_modelo_final2[(est_modelo_final2$op=="~"),]
    est_modelo_final_reg_aux$interp_pred<-NA
    est_modelo_final_reg_aux$interp_pred[is.na(est_modelo_final_reg_aux$interp_pred) & est_modelo_final_reg_aux$op=="~"]<-paste("O aumento de uma unidade na satisfacao da variavel latente",
                                                                                                                                est_modelo_final_reg_aux$label_rhs,"implica no aumento de",
                                                                                                                                round(est_modelo_final_reg_aux$est,4), "unidades na variavel observada",
                                                                                                                                est_modelo_final_reg_aux$label_lhs)
    
    est_modelo_final_reg_val<-rbind(est_modelo_final_val_aux,
                                    est_modelo_final_reg_aux)
    
    est_modelo_final_reg_val<-est_modelo_final_reg_val %>%
      select(lhs,op,rhs,est,pvalue,interp_pred)
    
    r2_bancofim<-est_modelo_final2[est_modelo_final2$op=="r2",]
    r2_bancofim<-r2_bancofim[,colSums(is.na(r2_bancofim))<nrow(r2_bancofim)]
    r2_bancofim_v1<-r2_bancofim[!duplicated(as.list(r2_bancofim))]
    
    
    
    est_modelo_final3<-est_modelo_final2[est_modelo_final2$op!="r2",]
    
    ##########
    # impacto
    
    
    est_modelo_impacto_aux<-est_modelo_final2[(est_modelo_final2$op=="~" | est_modelo_final2$op=="=~"),
                                              c("lhs", "label_lhs", "op", "rhs", "label_rhs","std.all")]
    
    est_modelo_impacto<-merge(est_modelo_impacto_aux,
                              df.top2box,
                              by.x="rhs",
                              by.y="variaveis", all.x=T)
    
    est_modelo_impacto<-est_modelo_impacto[c("lhs", "label_lhs", "op", "rhs", "label_rhs","std.all","T2B_result")]
    est_modelo_impacto$T2B_result<-as.numeric(as.character(est_modelo_impacto$T2B_result))
    # }
    #}
    ################################################################################################################################################################################################################
    ### Exportar para um excel
    ################################################################################################################################################################################################################
    
    # Create an Excel workbook. 
    # Both .xls and .xlsx file formats can be used.
    #setwd(export)
    
    #### vamos criar uma nova pasta das quebras
    
    #dir.create(nivel_label)
    #pasta_output<-paste(export, nivel_label,sep="\\")
    fileName<-paste("FGTeste",paste(nivel_label,"resultados_SEM.xlsx",sep="_"))
    #fileXls <- paste(pasta_output,fileName,sep='\\')
    fileXls <- fileName
    unlink(fileXls, recursive = FALSE, force = FALSE)
    exc <- loadWorkbook(fileXls, create = TRUE)
    
    ## Estimativas na planilha estimativas
    
    createSheet(exc,'Estimativas')
    saveWorkbook(exc)
    input <- est_modelo_final_reg_val
    writeWorksheet(exc, input, sheet ='Estimativas', startRow = 1, startCol = 2)
    saveWorkbook(exc)
    
    ### Analise Fatorial Confirmatoria Diagnostico Stepwise
    fileGraph <- paste("FGTeste",paste(nivel_label,'graph.png',sep="_"))
    png(filename = fileGraph, width = 800, height = 600)
    par(mfrow=c(2,3))
    plot(razao_chisq,ylim=c(0,70), main="Quanto mais perto de 5 melhor", type="l")
    abline(h=5,col="blue")
    plot(aic, main="Quanto menor, melhor", type="l")
    plot(rmsea, ylim=c(0.01,0.1),main="Quanto menor melhor", type="l")
    abline(h=0.085,col="red")
    plot(cfi, ylim=c(0.5,1),main="Quanto maior melhor", type="l")
    abline(h=0.7,col="red")
    plot(alfa_cronbach, main="Quanto maior melhor", 
         ylim=c(0.5,1), type="l")
    plot(ave_min, main="Quanto maior melhor", ylim=c(0.1,1), type="l")
    abline(h=0.5,col="red")
    invisible(dev.off())
    drop_upload(fileGraph, path = "FPB")
    
    
    ### Estatisticas resumo modelo
    resumo_modelo<-fitmeasures(modelo_final_final, c("rmsea", "cfi"))
    createSheet(exc,'EstatResumo')
    saveWorkbook(exc)
    input <- as.data.frame(cbind(round(resumo_modelo,4), 
                                 rownames(resumo_modelo)))
    input<-rbind(input,n=round(n_obs,1))
    
    input$estat<-rownames(input)
    writeWorksheet(exc, input, sheet ='EstatResumo',header=T,rownames = T, startRow = 1, startCol = 2)
    saveWorkbook(exc)
    
    ### R squared
    
    createSheet(exc,'Rsquared')
    saveWorkbook(exc)
    input <- r2_bancofim_v1
    writeWorksheet(exc, input, sheet ='Rsquared',header=T,rownames = T, startRow = 1, startCol = 2)
    saveWorkbook(exc)
    #drop_upload(fileName)
    #### Reliability CFA
    
    createSheet(exc,'ReliabilityCFA')
    saveWorkbook(exc)
    rel_aux<-reliability(modelo_partial[[melhor]])
    input <- data.frame(Estat=row.names(rel_aux),rel_aux) 
    writeWorksheet(exc, input, 
                   sheet ='ReliabilityCFA',header=T,rownames = T, startRow = 1, startCol = 2)
    
    
    
    
    input2<-data.frame(reliabilityL2(modelo_final_final,segunda_ordem))
    input3<-data.frame(cbind(CR=row.names(input2),segunda_ordem,est=(input2[,1])))
    input3$est<-as.numeric(as.character(input3$est))
    
    writeWorksheet(exc, input3, 
                   sheet ='ReliabilityCFA',header=T,rownames = T, startRow = nrow(input)+7, startCol = 2)
    
    modelo_partial = NULL
    ##############
    saveWorkbook(exc)
    
    
    
    ### Impactos
    
    createSheet(exc,'Impactos')
    saveWorkbook(exc)
    input <- est_modelo_impacto[order(est_modelo_impacto$lhs),]
    
    
    writeWorksheet(exc, input, sheet ='Impactos',header=T,rownames = T, startRow = 1, startCol = 2)
    saveWorkbook(exc)
    
    ## Diferenca entre variaveis entrada/saida
    difference<-NULL
    difference <- getDifference(difference)
    difference2<-NULL
    # difference2 <-data.frame(codigo=unique(unlist(list(difference, variaveis_PROBLEMA_MATRIZ))))
    # difference2$codigo<-as.character(difference2$codigo)
    difference2<-as.data.frame(difference)
    difference2$difference<-as.character(difference2$difference)
    var_labels2<-var_labels
    var_labels2$codigo<-as.character(var_labels2$codigo)
    
    
    var_fora<-merge(difference2, var_labels2, all.x=T, by.x="difference", by.y = "codigo")
    colnames(var_fora)[1]<-("VariaveisDeFora")
    
    createSheet(exc,'VarFora')
    saveWorkbook(exc)
    input <- var_fora
    writeWorksheet(exc, input, sheet ='VarFora',header=T,rownames = T, startRow = 1, startCol = 2)
    saveWorkbook(exc)
    drop_upload(fileName, path = "FPB")
    
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

save.image("FG.RData")
history("FG.Rhistory")

Fim_Sintaxe_SEM<-Sys.time()
tempo_total<-Fim_Sintaxe_SEM - TempoTotal
Fim_Sintaxe_SEM - TempoTotal
tempos2 <- imprimirTempoPonderacao(tempo_total, "Tempo total: ")
indiceTempos2 <- indiceTempos2 +1
################################################################################################################################################################################################################
################################################################################################################################################################################################################
################################################################################################################################################################################################################
################################################################################################################################################################################################################
################################################################################################################################################################################################################
################################################################################################################################################################################################################