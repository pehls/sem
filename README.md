# sem - Projeto de Complex Survey Structural Equation Modeling 
## Roda em uma VM ubuntu 17.10 no google cloud platform

# GCP Config 

## Machine type 
#### n1-highmem-4 (4 vCPUs, 26 GB de memória),
#### hd 10 gb
## Debian Config on GCP
#### sudo apt-get install r-base
#### sudo apt-get install openjdk-8-jdk
#### sudo apt-get install r-cran-rjava
#### sudo apt-get install curl
#### sudo apt-get install libssl-dev
#### sudo apt-get install libcurl4-openssl-dev
#### sudo apt-get install libxml2-dev
#### sudo apt-get install libgsl0ldbl
#### sudo apt-get install gsl-bin libgsl0-de
## To open R, simply write "R"
## if process was killed sometime, open  /var/log/kern.log to examinate it;
## or use dmesg to examinate the kernel log, too.
