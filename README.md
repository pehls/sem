# sem - Projeto de Complex Survey Structural Equation Modeling

Objeto de Estudo do meu TCC em Sistemas de Informação, na UniRitter.

Visa analisar o impacto da remoção de outliers e da utilização de Bootstrap nas Equações Estruturais, no que tange a métricas e também na velocidade e requerimentos de hardware para a execução do mesmo.

## Originalmente, roda em uma VM Debian 9 no google cloud platform

# GCP Config

## Machine type
#### n1-highmem-8 (8 vCPUs, 52 GB de memória),
#### hd 10 gb
## Debian Config on GCP
```
 sudo apt-get install r-base
 sudo apt-get install openjdk-8-jdk
 sudo apt-get install r-cran-rjava
 sudo apt-get install curl
 sudo apt-get install libssl-dev
 sudo apt-get install libcurl4-openssl-dev
 sudo apt-get install libxml2-dev
 sudo apt-get install libgsl0ldbl
 sudo apt-get install gsl-bin libgsl0-de
 sudo nano /etc/apt/sources.list
```
### Dentro do arquivo adicionar (ou seja, editando o mesmo):
```
 deb http://cran.rstudio.com/bin/linux/debian stretch-cran35/
```
 ctrl o enter (para salvar)
 ctrl x (para sair)
 ```
 sudo apt-get update (provavelmente dará problema com assinatura do pacote, não tem problema :) )
 sudo apt-get install r-base
 ```
#### Para abrir o R, apenas digite "R"
#### Se o processo for 'killed', abra o  /var/log/kern.log para examinar, ou use dmesg para tal.
