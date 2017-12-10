# sem - Projeto de Complex Survey Structural Equation Modeling 
## Roda em uma VM ubuntu 17.10 no google cloud platform

# GCP Config 

## Machine type 
\n n1-highmem-8 (8 vCPUs (analyzing the needs/quantities), 52 GB de memória),
\n ssd 50 gb
## CPU Platform
\n Intel Broadwell
## Zone
\n southamerica-east1-c
## Debian Config on GCP
\n sudo apt-get install r-base
\n sudo apt-get install openjdk-8-jdk
\n sudo apt-get install openjdk-8-source #opcional
\n sudo apt-get install r-cran-rjava
\n sudo apt-get install curl
\n sudo apt-get install libssl-dev
\n sudo apt-get install libcurl4-openssl-dev
\n sudo apt-get install libxml2-dev
\n sudo apt-get install libgsl0ldbl
\n sudo apt-get install gsl-bin libgsl0-de
## To open R, simply write "R"