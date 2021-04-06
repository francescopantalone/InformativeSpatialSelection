dir1<-try(file.path(Mydirectories::googledrive.directory(),"Travail/Recherche/Travaux","InformativeSpatialSelection"))
dir2<-try(file.path("/Users/francesco/Google Drive/InformativeSpatialSelection"))
if(dir.exists(dir1)){setwd(dir1)}
if(dir.exists(dir2)){setwd(dir2)}
library(InformativeSpatialSelection)
X<-sapply(list.files("R",full.names = T),function(x){try(source(x))})
which(is.element(sapply(X,class),"try-error"))
#Simulation of data
##creation of 100 replicates of Y and epsilon
if(FALSE){
  set.seed(1)
  nrep=3000
  .scale=10
  x1<-generatex_i(201)
  Y3000<-generateYonsquaregridU(x1=x1,x2=x1, mu=0,nrep=nrep)
  epsilon3000<-generateYonsquaregridU(x1=x1,x2=x1, mu=0,nrep=nrep)
  save(Y3000, 
       file = file.path(getwd(), "data", "Y3000.rda"))
  #load(file.path(getwd(), "data", "Y3000.rda"))
  save(epsilon3000, 
       file = file.path(getwd(), "data", "epsilon3000.rda"))
  #load(file.path(getwd(), "data", "epsilon3000.rda"))
}


##create datasets figure<x>_data <x>=1...6
figure1_data<-Simulate_figure1_data();
save(figure1_data, 
     file = file.path(getwd(), "data", "figure1_data.rda"))

figure1_data<-get(load(file.path(getwd(), "data", "figure1_data.rda")))
figure2_data<-Simulate_figure2_data(figure1_data=figure1_data);
save(figure2_data, 
     file = file.path(getwd(), "data", "figure2_data.rda"))

figure1_data<-get(load(file.path(getwd(), "data", "figure1_data.rda")))
figure3_data<-Simulate_figure3_data(figure1_data=figure1_data);
save(figure3_data, file = file.path(getwd(), "data", "figure3_data.rda"))

figure1_data<-get(load(file.path(getwd(), "data", "figure1_data.rda")))
figure4_data<-Simulate_figure4_data(figure1_data=figure1_data)
save(figure4_data, file = file.path(getwd(), "data", "figure4_data.rda"))

figure1_data<-get(load(file.path(getwd(), "data", "figure1_data.rda")))
figure6_data<-Simulate_figure6_data(figure1_data=figure1_data)
save(figure6_data, file = file.path(getwd(), "data", "figure6_data.rda"))

# Creation of figures
#source('R/3.1.figure1.R')
figure1(figure1_data=eval(get(load(file.path(getwd(), "data", "figure1_data.rda")))))
figure1(
  pdfexportfolder=file.path(getwd(),"latex/fig"),
        pngconvert = TRUE,
        figure1_data=eval(get(load(file.path(getwd(), "data", "figure1_data.rda")))))
fs::file_show(file.path(getwd(),"latex/fig/figure1.pdf"))
fs::file_show(file.path(getwd(),"latex/fig/figure1.png"))

source('R/3.2.figure2.R')
figure2(figure2_data=eval(get(load(file.path(getwd(), "data", "figure2_data.rda")))))
figure2(
  pdfexportfolder=file.path(getwd(),"latex/fig"),
  pngconvert = TRUE,
  figure2_data=eval(get(load(file.path(getwd(), "data", "figure2_data.rda")))))
fs::file_show(file.path(getwd(),"latex/fig/figure2.pdf"))
#fs::file_show(file.path(getwd(),"latex/fig/figure2.png"))

source('R/3.3.figure3.R')
figure3(
  pdfexportfolder=file.path(getwd(),"latex/fig"),
  pngconvert = TRUE)
fs::file_show(file.path(getwd(),"latex/fig/figure3.pdf"))
#fs::file_show(file.path(getwd(),"latex/fig/figure3.png"))

source('R/figure4.R')
figure4(
  pdfexportfolder=file.path(getwd(),"latex/fig"),
  pngconvert = TRUE)
fs::file_show(file.path(getwd(),"latex/fig/figure4.pdf"))
fs::file_show(file.path(getwd(),"latex/fig/figure4.png"))

figure5(
  pdfexportfolder=file.path(getwd(),"latex/fig"),
  pngconvert = TRUE)
fs::file_show(file.path(getwd(),"latex/fig/figure5.pdf"))
fs::file_show(file.path(getwd(),"latex/fig/figure5.png"))

