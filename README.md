# sem - Projeto de Complex Survey Structural Equation Modeling

-Objeto de Estudo do meu TCC em Sistemas de Informação, na UniRitter.

-Visa analisar o impacto da remoção de outliers e da utilização de Bootstrap nas Equações Estruturais, no que tange a métricas e também na velocidade e requerimentos de hardware para a execução do mesmo;

-Através de Métricas comuns às Equações Estruturais, como o CFI, RMSEA e um teste de Diferenças significativas dos coeficientes, via Teste Z, implementado via campos calculados no Tableau Desktop, verificamos que a utilização de um bootstrap com o mínimo de 400 linhas aleatórias da nossa base total não modifica de forma significativa os mesmos, assim como mantém o CFI e a raiz do erro médio quadrático em um nível adequado segundo a literatura (Hair, Black, Babin, Anderson @ Multivariate Data Analysis, 2014);

-Identificada a possibilidade de diminuir o tempo de execução para entre 8 e 11 minutos (por vez), diminuindo a memória RAM alocada para o R de um pico de 12GB para cerca de 7GB, tornando possível a utilização de um computador local ao invés de uma máquina virtual, conforme o projeto original;

-Necessita de maiores testes em mercados diferentes, devido á proximidade dos mercados das IES utilizadas, bem como testes com tamanhos de bases diferentes, para comprovar os resultados.

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
