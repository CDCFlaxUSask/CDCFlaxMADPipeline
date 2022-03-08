#The installs for the code see below
#install.packages("agricolae",destdir = "d:/download/RDownload")
#install.packages("Rtools",destdir = "d:/download/RDownload")
#install.packages("asreml",destdir = "d:/download/RDownload")
#install.packages("readxl",destdir = "d:/download/RDownload")
#install.packages("tidyverse",destdir = "d:/download/RDownload")
#install.packages("data.table",destdir = "d:/download/RDownload")
#install.packages("ggplot2",destdir = "d:/download/RDownload")
#install.packages("jsonlite",destdir = "d:/download/RDownload")
#install.packages("stringr",destdir = "d:/download/RDownload")
#install.packages("D:/RStuff/ASReml/asreml_4.1.0.106.zip", repos = NULL, type = "win.binary")
#install.packages("rlist")
#asreml.license.activate()#Please enter your activation code (RET or 0 to exit):XXXX-XXXX-XXXX-XXXX
#  install.packages("sasLM")

#install.packages("knitr") # Run if you don't have this one
#install.packages("sasLM", repos="http://r.acr.kr")
#install.packages("xlsx",destdir = "d:/download/RDownload")

########################################################################################################

library(stringr)
library(asreml)
library("agricolae") ##not actually necessary for LSD 5% calcs but can make things easier if the data is balanced. 
library("readxl") ## just to read the data
library("tidyverse") ## just to read the data
library(stringr)
library(rlist)

library(dplyr)
library(tidyr)

options(max.print = 999999)
options(tibble.print_max=50)

#this opens up a file upload prompt
DataFilename<-file.choose(new = FALSE)
DataFilenameConverted<-str_replace(DataFilename, ".txt", "_Converted.csv")
DataFilenameConvertedNew<-str_replace(DataFilename, ".txt", "_Converted_New.csv")
DataFileDir<-str_replace(DataFilename, "BxN_Preston_Wilt_Data_2019_for_MAD.txt", "")
DataFileDircsv<-paste(str_replace(DataFilename, "BxN_Preston_Wilt_Data_2019_for_MAD.txt", ""),"csv\\",sep="")

#get the OS we are on windows ???
OS<-Sys.info()["sysname"]

#if windows get the default dir otherwise its UN*X/MAC
if (.Platform$OS.type == "windows"){
  dirend<-sapply(gregexpr("\\\\", DataFilename), tail, 1)
}else {
  dirend<-sapply(gregexpr("\\/", DataFilename), tail, 1)
}
diris<-substr(DataFilename, 1, dirend)

#CSV DIR path
csvdir<-paste(diris,"txt\\", sep = "")
txtdir<-paste(diris,"csv\\", sep = "")
workdir<-paste(diris,"work\\", sep = "")

#if the CSV subdir s\doesnt exist create it
if (file.exists(csvdir)== FALSE)
{
  dir.create(file.path(csvdir))
}
#if the WORK subdir s\doesnt exist create it
if (file.exists(workdir)== FALSE)
{
  dir.create(file.path(workdir))
}

#set working directory
setwd(diris)

#Read in the data from the selected file
#data<-read_excel(DataFilename)
data<-read_tsv(DataFilename)
               



#Record	Plot	Row	Column	Cp	Csp	Entry	Year	Location	Genotype
header<-names( data )#names(data)

correct<-grepl('Record, Plot, Row, Column, Cp, Csp, Entry, Year, Location, Genotype', toString(header))

if (correct) 
  {
  
  DataConverted <-tidyr::pivot_longer(data, STD:W3, names_to = "Trait", values_to = "Value")
  
  Locations<-unique(DataConverted$Location)
  Years<-unique(DataConverted$Year)
  Traits<-unique(DataConverted$Trait)
  
  
  write.csv(DataConverted, file = DataFilenameConverted,row.names=TRUE)
  data2<-read_csv(DataFilenameConverted)
  
  data2$Whole_plot <- "."

  
  for( currentRow in 1:nrow(data2))
    {
      if ( data2[currentRow, 6] == 1 || data2[currentRow, 7] > 0) 
      {
        data2[currentRow, ncol(data2)]  <- paste(data2[currentRow,4],"+",data2[currentRow,5])
      }
    }
  
  
  Locations<-unique(data2$Location)
  Years<-unique(data2$Year)
  Traits<-unique(data2$Trait)

  
  write.csv(data2, file = DataFilenameConvertedNew,row.names=TRUE)
  
  for( Location in Locations)
  {
   for( Year in Years)
   {
    print(paste(Year,Location))

  
  
    mad_phenotypic_data <- data2
    MDA_subplot_controls  <- mad_phenotypic_data %>% filter(Csp > 0) 
    
    
    
    MDA_plot_subplot_controls  <- mad_phenotypic_data %>% inner_join(MDA_subplot_controls,  by = c('Whole_plot'='Whole_plot', 'Year'='Year', 'Location'='Location', 'Trait'='Trait')) %>% filter(Whole_plot != ".") 
    MDA_plot_subplot_controls  <- MDA_plot_subplot_controls[1:14] %>% distinct()
    MDA_plot_subplot_controls <- MDA_plot_subplot_controls[order(MDA_plot_subplot_controls$Year, MDA_plot_subplot_controls$Location, MDA_plot_subplot_controls$Trait, MDA_plot_subplot_controls$Whole_plot),]
  
    MDA_plot_controls <- mad_phenotypic_data   %>% filter(Cp == 1) %>% filter(Csp == 0) #%>% filter(Trait == 'V1') 
 
    write.csv(mad_phenotypic_data, file = paste(DataFileDir,"mad_phenotypic_data",".csv"),row.names=TRUE)
    write.csv(MDA_plot_subplot_controls, file = paste(DataFileDir,"MDA_plot_subplot_controls",".csv"),row.names=TRUE)
    write.csv(MDA_plot_controls, file = paste(DataFileDir,"MDA_plot_controls",".csv"),row.names=TRUE)
  
    require(sasLM)

    d1 = read.csv(" MDA_plot_controls .csv")
    d1 = af(d1, c("Year", "Location", "Trait")) 
    f1 = Value ~ factor(Row) + factor(Column)
    x4<-BY(GLM, f1, d1, By="Trait")
  
    d2 = read.csv(" MDA_plot_subplot_controls .csv")
    d2 = af(d2, c("Year", "Location", "Trait")) 
    f2 = Value.x ~ factor(Whole_plot) + factor(Genotype.x)
    x5<-BY(GLM, f2, d2, By="Trait")
  
  
    #capture.output(x4$STD$ANOVA, file = paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".txt",sep=""))
    #capture.output(x4$STD$ANOVA[2][1], file = paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".txt",sep=""),append = TRUE)
    #HeaderTXT <- paste("Year,Location,Type,Df,Sum Sq,Mean Sq,F value,Pr(>F)",sep=",")
    #capture.output(HeaderTXT,file = paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".txt",sep=""))
    #stderr<-paste(Year,Location,"STD",x4$STD$ANOVA[[2]][[1]],x4$STD$ANOVA[[5]][[1]],x4$STD$ANOVA[[8]][[1]],x4$STD$ANOVA[[11]][[1]],x4$STD$ANOVA[[14]][[1]],sep=",")
    #capture.output(stderr,file = paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".txt",sep=""),append = TRUE)

    HeaderTXT <- paste("Year,Location,Type,Parameter,Df,Sum Sq,Mean Sq,F value,Pr(>F)",sep=",")
    
    filename2<-paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".csv",sep="")
    write(c("Year,Location,Type,Parameter,Df,Sum Sq,Mean Sq,F value,Pr(>F)"), file=filename2,append=FALSE)
    
    outline<-paste(Year,Location,"STD,ERROR_RXC",x4$STD$ANOVA[[2]][[1]],x4$STD$ANOVA[[5]][[1]],x4$STD$ANOVA[[8]][[1]],x4$STD$ANOVA[[11]][[1]],x4$STD$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"STD,ROW",x4$STD$`Type III`[[1]][[1]],x4$STD$`Type III`[[3]][[1]],x4$STD$`Type III`[[5]][[1]],x4$STD$`Type III`[[7]][[1]],x4$STD$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"STD,COL",x4$STD$`Type III`[[2]][[1]],x4$STD$`Type III`[[4]][[1]],x4$STD$`Type III`[[6]][[1]],x4$STD$`Type III`[[8]][[1]],x4$STD$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    

    filename2<-paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"V1,ERROR_RXC",x4$V1$ANOVA[[2]][[1]],x4$V1$ANOVA[[5]][[1]],x4$V1$ANOVA[[8]][[1]],x4$V1$ANOVA[[11]][[1]],x4$V1$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"V1,ROW",x4$V1$`Type III`[[1]][[1]],x4$V1$`Type III`[[3]][[1]],x4$V1$`Type III`[[5]][[1]],x4$V1$`Type III`[[7]][[1]],x4$V1$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"V1,COL",x4$V1$`Type III`[[2]][[1]],x4$V1$`Type III`[[4]][[1]],x4$V1$`Type III`[[6]][[1]],x4$V1$`Type III`[[8]][[1]],x4$V1$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"V2,ERROR_RXC",x4$V2$ANOVA[[2]][[1]],x4$V2$ANOVA[[5]][[1]],x4$V2$ANOVA[[8]][[1]],x4$V2$ANOVA[[11]][[1]],x4$V2$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"V2,ROW",x4$V2$`Type III`[[1]][[1]],x4$V2$`Type III`[[3]][[1]],x4$V2$`Type III`[[5]][[1]],x4$V2$`Type III`[[7]][[1]],x4$V2$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"V2,COL",x4$V2$`Type III`[[2]][[1]],x4$V2$`Type III`[[4]][[1]],x4$V2$`Type III`[[6]][[1]],x4$V2$`Type III`[[8]][[1]],x4$V2$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"V3,ERROR_RXC",x4$V3$ANOVA[[2]][[1]],x4$V3$ANOVA[[5]][[1]],x4$V3$ANOVA[[8]][[1]],x4$V3$ANOVA[[11]][[1]],x4$V3$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"V3,ROW",x4$V3$`Type III`[[1]][[1]],x4$V3$`Type III`[[3]][[1]],x4$V3$`Type III`[[5]][[1]],x4$V3$`Type III`[[7]][[1]],x4$V3$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"V3,COL",x4$V3$`Type III`[[2]][[1]],x4$V3$`Type III`[[4]][[1]],x4$V3$`Type III`[[6]][[1]],x4$V3$`Type III`[[8]][[1]],x4$V3$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"W1,ERROR_RXC",x4$W1$ANOVA[[2]][[1]],x4$W1$ANOVA[[5]][[1]],x4$W1$ANOVA[[8]][[1]],x4$W1$ANOVA[[11]][[1]],x4$W1$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"W1,ROW",x4$W1$`Type III`[[1]][[1]],x4$W1$`Type III`[[3]][[1]],x4$W1$`Type III`[[5]][[1]],x4$W1$`Type III`[[7]][[1]],x4$W1$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"W1,COL",x4$W1$`Type III`[[2]][[1]],x4$W1$`Type III`[[4]][[1]],x4$W1$`Type III`[[6]][[1]],x4$W1$`Type III`[[8]][[1]],x4$W1$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"W2,ERROR_RXC",x4$W2$ANOVA[[2]][[1]],x4$W2$ANOVA[[5]][[1]],x4$W2$ANOVA[[8]][[1]],x4$W2$ANOVA[[11]][[1]],x4$W2$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"W2,ROW",x4$W2$`Type III`[[1]][[1]],x4$W2$`Type III`[[3]][[1]],x4$W2$`Type III`[[5]][[1]],x4$W2$`Type III`[[7]][[1]],x4$W2$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"W2,COL",x4$W2$`Type III`[[2]][[1]],x4$W2$`Type III`[[4]][[1]],x4$W2$`Type III`[[6]][[1]],x4$W2$`Type III`[[8]][[1]],x4$W2$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"W3,ERROR_RXC",x4$W3$ANOVA[[2]][[1]],x4$W3$ANOVA[[5]][[1]],x4$W3$ANOVA[[8]][[1]],x4$W3$ANOVA[[11]][[1]],x4$W3$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"W3,ROW",x4$W3$`Type III`[[1]][[1]],x4$W3$`Type III`[[3]][[1]],x4$W3$`Type III`[[5]][[1]],x4$W3$`Type III`[[7]][[1]],x4$W3$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"W3,COL",x4$W3$`Type III`[[2]][[1]],x4$W3$`Type III`[[4]][[1]],x4$W3$`Type III`[[6]][[1]],x4$W3$`Type III`[[8]][[1]],x4$W3$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    
    HeaderTXT <- paste("Year,Location,Type,Parameter,Df,Sum Sq,Mean Sq,F value,Pr(>F)",sep=",")
    
    filename2<-paste(csvdir,"unadjusted_subplot_anova_stats",Year,Location,".csv",sep="")
    subplot_anova_file<-filename2
    
    write(c("Year,Location,Type,Parameter,Df,Sum Sq,Mean Sq,F value,Pr(>F)"), file=filename2,append=FALSE)
    
    outline<-paste(Year,Location,"STD,ERROR_WG",x5$STD$ANOVA[[2]][[1]],x5$STD$ANOVA[[5]][[1]],x5$STD$ANOVA[[8]][[1]],x5$STD$ANOVA[[11]][[1]],x5$STD$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"STD,Whole_Plot",x5$STD$`Type III`[[1]][[1]],x5$STD$`Type III`[[3]][[1]],x5$STD$`Type III`[[5]][[1]],x5$STD$`Type III`[[7]][[1]],x5$STD$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"STD,Genotype",x5$STD$`Type III`[[2]][[1]],x5$STD$`Type III`[[4]][[1]],x5$STD$`Type III`[[6]][[1]],x5$STD$`Type III`[[8]][[1]],x5$STD$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    
    filename2<-paste(csvdir,"unadjusted_subplot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"V1,ERROR_WG",x5$V1$ANOVA[[2]][[1]],x5$V1$ANOVA[[5]][[1]],x5$V1$ANOVA[[8]][[1]],x5$V1$ANOVA[[11]][[1]],x5$V1$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"V1,Whole_Plot",x5$V1$`Type III`[[1]][[1]],x5$V1$`Type III`[[3]][[1]],x5$V1$`Type III`[[5]][[1]],x5$V1$`Type III`[[7]][[1]],x5$V1$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"V1,Genotype",x5$V1$`Type III`[[2]][[1]],x5$V1$`Type III`[[4]][[1]],x5$V1$`Type III`[[6]][[1]],x5$V1$`Type III`[[8]][[1]],x5$V1$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_subplot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"V2,ERROR_WG",x5$V2$ANOVA[[2]][[1]],x5$V2$ANOVA[[5]][[1]],x5$V2$ANOVA[[8]][[1]],x5$V2$ANOVA[[11]][[1]],x5$V2$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"V2,Whole_Plot",x5$V2$`Type III`[[1]][[1]],x5$V2$`Type III`[[3]][[1]],x5$V2$`Type III`[[5]][[1]],x5$V2$`Type III`[[7]][[1]],x5$V2$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"V2,Genotype",x5$V2$`Type III`[[2]][[1]],x5$V2$`Type III`[[4]][[1]],x5$V2$`Type III`[[6]][[1]],x5$V2$`Type III`[[8]][[1]],x5$V2$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_subplot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"V3,ERROR_WG",x5$V3$ANOVA[[2]][[1]],x5$V3$ANOVA[[5]][[1]],x5$V3$ANOVA[[8]][[1]],x5$V3$ANOVA[[11]][[1]],x5$V3$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"V3,Whole_Plot",x5$V3$`Type III`[[1]][[1]],x5$V3$`Type III`[[3]][[1]],x5$V3$`Type III`[[5]][[1]],x5$V3$`Type III`[[7]][[1]],x5$V3$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"V3,Genotype",x5$V3$`Type III`[[2]][[1]],x5$V3$`Type III`[[4]][[1]],x5$V3$`Type III`[[6]][[1]],x5$V3$`Type III`[[8]][[1]],x5$V3$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_subplot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"W1,ERROR_WG",x5$W1$ANOVA[[2]][[1]],x5$W1$ANOVA[[5]][[1]],x5$W1$ANOVA[[8]][[1]],x5$W1$ANOVA[[11]][[1]],x5$W1$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"W1,Whole_Plot",x5$W1$`Type III`[[1]][[1]],x5$W1$`Type III`[[3]][[1]],x5$W1$`Type III`[[5]][[1]],x5$W1$`Type III`[[7]][[1]],x5$W1$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"W1,Genotype",x5$W1$`Type III`[[2]][[1]],x5$W1$`Type III`[[4]][[1]],x5$W1$`Type III`[[6]][[1]],x5$W1$`Type III`[[8]][[1]],x5$W1$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_subplot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"W2,ERROR_WG",x5$W2$ANOVA[[2]][[1]],x5$W2$ANOVA[[5]][[1]],x5$W2$ANOVA[[8]][[1]],x5$W2$ANOVA[[11]][[1]],x5$W2$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"W2,Whole_Plot",x5$W2$`Type III`[[1]][[1]],x5$W2$`Type III`[[3]][[1]],x5$W2$`Type III`[[5]][[1]],x5$W2$`Type III`[[7]][[1]],x5$W2$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"W2,Genotype",x5$W2$`Type III`[[2]][[1]],x5$W2$`Type III`[[4]][[1]],x5$W2$`Type III`[[6]][[1]],x5$W2$`Type III`[[8]][[1]],x5$W2$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)
    
    filename2<-paste(csvdir,"unadjusted_subplot_anova_stats",Year,Location,".csv",sep="")
    outline<-paste(Year,Location,"W3,ERROR_WG",x5$W3$ANOVA[[2]][[1]],x5$W3$ANOVA[[5]][[1]],x5$W3$ANOVA[[8]][[1]],x5$W3$ANOVA[[11]][[1]],x5$W3$ANOVA[[14]][[1]],sep=",")
    outline1<-paste(Year,Location,"W3,Whole_Plot",x5$W3$`Type III`[[1]][[1]],x5$W3$`Type III`[[3]][[1]],x5$W3$`Type III`[[5]][[1]],x5$W3$`Type III`[[7]][[1]],x5$W3$`Type III`[[9]][[1]],sep=",")
    outline2<-paste(Year,Location,"W3,Genotype",x5$W3$`Type III`[[2]][[1]],x5$W3$`Type III`[[4]][[1]],x5$W3$`Type III`[[6]][[1]],x5$W3$`Type III`[[8]][[1]],x5$W3$`Type III`[[10]][[1]],sep=",")
    write(outline, file=filename2,append=TRUE)
    write(outline1, file=filename2,append=TRUE)
    write(outline2, file=filename2,append=TRUE)

   }
  }
  

  datap1 <- read.csv(DataFilenameConvertedNew)[ ,3:14]  
  DataFilenameConvertedNewOut <- str_replace(DataFilenameConvertedNew, ".csv", "_Out.csv")
  datap1 <- datap1   %>% filter(Trait != '.') 
  datap1$ID <- paste(datap1$Year ,datap1$Location,datap1$Genotype, sep="_")#create an environment attribute in the dataframe
#  datap1$IDROW <- paste(datap1$Year ,datap1$Location,datap1$Genotype,datap1$Row, sep="_")#create an environment attribute in the dataframe
#  datap1$IDCOL <- paste(datap1$Year ,datap1$Location,datap1$Genotype,datap1$Col, sep="_")#create an environment attribute in the dataframe
  datap1$IDROW <- paste(datap1$Year ,datap1$Location,datap1$Trait,datap1$Row, sep="_")#create an environment attribute in the dataframe
  datap1$IDCOL <- paste(datap1$Year ,datap1$Location,datap1$Trait,datap1$Col, sep="_")#create an environment attribute in the dataframe
  datap1$IDX <- paste(datap1$Year ,datap1$Location,datap1$Trait, sep="_")#create an environment attribute in the dataframe
#  datap1$IDROWCOL <- paste(datap1$Year ,datap1$Location,datap1$Genotype,datap1$Row,datap1$Col, sep="_")#create an environment attribute in the dataframe
  datap1$IDROWCOL <- paste(datap1$Year ,datap1$Location,datap1$Trait,datap1$Row,datap1$Col, sep="_")#create an environment attribute in the dataframe
  
    DataConverted$CTRLGENOTYPE <- 0#create an environment attribute in the dataframe
  
  
  if (DataConverted$Cp == 1 || DataConverted$Csp > 0) {
    DataConverted$CTRLGENOTYPE <- 1;
  }
  
  write.csv(DataConverted, file = paste(DataFileDir,"DataConverted",".csv"),row.names=FALSE)
  write.csv(DataConverted, file = paste(DataFileDircsv,"aaaaDataConverted",".csv"),row.names=FALSE)
  
  
  write.csv(data, file = paste(DataFileDir,"data",".csv"),row.names=FALSE)
  
  head(datap1)
  write.csv(datap1, file = DataFilenameConvertedNewOut,row.names=FALSE)  
  cp1 <- datap1 %>% filter(Cp == "1")
  ncp1 <- datap1 %>% filter(Cp != "1")
  
  write.csv(cp1, file = paste(DataFileDir,"cp1s",".csv"),row.names=FALSE)
  write.csv(ncp1, file = paste(DataFileDir,"ncp1",".csv"),row.names=FALSE)
  
  
  #-------------output
  p_row_means2 <- aggregate( Value ~ IDROW , cp1, mean )#row_means
  p_col_means2 <- aggregate( Value ~ IDCOL , cp1, mean )#col_means
  #m1_row_col_values <-paste(cp1$year,"-",cp1$Location,"-",cp1$Trait),cp1$row,cp1$Column,cp1$Value#row_col_values
  #p_row_col_values2 <- cp1 %>%  unite(ID, "Year", "Location", "Trait", sep = "-") %>%  select(ID, Row, Column, Value)
  p_row_col_values2 <- cp1 %>%  unite(IDROWCOL, "Year", "Location", "Trait", "Row", "Column", sep = "_") %>%  select(IDROWCOL, Value)
  p_row_col_testplot_means2 <- aggregate( Value ~ IDROWCOL , ncp1, mean )#row_col_testplot_means
  m1_means <- aggregate( Value ~ ID , cp1, mean )#total_means
 p_total_means <- aggregate( Value ~ IDX , cp1, mean )#total_means
  row_ids <- unique(datap1$Row)#datap1$IDROW#row_ids
  col_ids <- unique(datap1$Column)#datap1$IDCOL#col_ids
  ControlGenoType <- DataConverted$CTRLGENOTYPE#control_genotypes
  datap1cg <- filter(data, (data$Csp > 0) | (data$Cp == 1))
  

    
  head(cp1)
  head(datap1)
  
  unique(cp1$Genotype)
  
  #debug
  write.csv(datap1, file = paste(DataFileDir,"datap1",".csv"),row.names=FALSE)
  #debug
  
  p_row_col_values <- p_row_col_values2 
  #p_row_col_values <- p_row_col_values2[order(p_row_col_values$IDROWCOL),] #order(foo$V1)
  p_row_col_testplot_means <- p_row_col_testplot_means2
  #p_row_col_testplot_means <- p_row_col_testplot_means2[order(p_row_col_testplot_means$IDROWCOL),]
  
  p_row_col_values[order(p_row_col_values$IDROWCOL),]
  p_row_col_testplot_means[order(p_row_col_testplot_means$IDROWCOL),]
  

  names(p_row_col_testplot_means)[names(p_row_col_testplot_means)=="Value"] <- "MeanValue"
  
  
  cp1$IDX <- paste(cp1$Year ,cp1$Location,cp1$Trait, sep="_")#create an environment attribute in the dataframe
  
  
  
  IDmatrix <- cp1[, c("IDROWCOL", "IDX")]
  ListIDS <- unique(IDmatrix$IDX)
  
  RegressionData <- merge( p_row_col_testplot_means, p_row_col_values, on='IDROWCOL')
  RegressionData <- merge( RegressionData, IDmatrix, on='IDROWCOL')
  
  rcoeffs <- data.frame()
  
  for (ListID in ListIDS) { #ListID <-"2019_Preston_V3"
    subd <- RegressionData[RegressionData$IDX == ListID, ]
    #print(ListID)
    #lm(Value ~ MeanValue, data = subd)
    res1<-lm(MeanValue ~ Value, data = subd)
    print(res1)
    summary(res1)
    op<-res1[["coefficients"]][["Value"]]
    newdf <- data.frame(ListID,op)
    # <- rbind(rcoeffs, newdf) 
    rcoeffs <- rbind(rcoeffs, newdf) 
  }
  write.csv(rcoeffs, file = paste(DataFileDircsv,"rcoeffs",".csv"),row.names=FALSE)
  
  unadjusted_plot_anova_stats    <- read.csv(paste(csvdir,"unadjusted_plot_anova_stats",Year,Location,".csv",sep="")) 
  unadjusted_plot_anova_stats[is.na(unadjusted_plot_anova_stats)] <- 99
  unadjusted_plot_anova_stats$Significance<-"ns"
  unadjusted_plot_anova_stats$errorrxc_ms<-0
  unadjusted_plot_anova_stats$rxc_f<-0
  unadjusted_plot_anova_stats$rxc_significance<-"ns"
  unadjusted_plot_anova_stats$fprob<-""
  unadjusted_plot_anova_stats$method<-""

  unadjusted_subplot_anova_stats <- read.csv(paste(csvdir,"unadjusted_subplot_anova_stats",Year,Location,".csv",sep=""))
  unadjusted_subplot_anova_stats[is.na(unadjusted_subplot_anova_stats)] <- 99
  unadjusted_subplot_anova_stats$Significance="ns"
  unadjusted_subplot_anova_stats$errorrxc_ms=0 
  unadjusted_subplot_anova_stats$rxc_f=0 
  unadjusted_subplot_anova_stats$rxc_significance="ns"
  unadjusted_subplot_anova_stats$fprob=""
  unadjusted_subplot_anova_stats$method=""
  
  ##calculate MS column for subplot
  
  for(i in 1:nrow(unadjusted_subplot_anova_stats)) { #i=12
    row
    if(unadjusted_subplot_anova_stats[i,9] != 99 && unadjusted_subplot_anova_stats[i,9] <= 0.05)
    {
      unadjusted_subplot_anova_stats[i,10]<-"*"
    }
    if(unadjusted_subplot_anova_stats[i,9] != 00 && unadjusted_subplot_anova_stats[i,9] <= 0.01)
    {
      unadjusted_subplot_anova_stats[i,10]<-"**"
    }
    
    if(unadjusted_plot_anova_stats[i,9] != 99 && unadjusted_plot_anova_stats[i,9] <= 0.05)
    {
      unadjusted_plot_anova_stats[i,10]<-"*"
    }
    if(unadjusted_plot_anova_stats[i,9] != 99 && unadjusted_plot_anova_stats[i,9] <= 0.01)
    {
      unadjusted_plot_anova_stats[i,10]<-"**"
    }
    
    unadjusted_subplot_anova_stats[i,11]<-unadjusted_subplot_anova_stats[i,6]/unadjusted_subplot_anova_stats[i,5]
    
    unadjusted_plot_anova_stats[i,11]<-unadjusted_plot_anova_stats[i,6]/unadjusted_plot_anova_stats[i,5]
    
    
    unadjusted_subplot_anova_stats[i,12] <- unadjusted_plot_anova_stats[i,11]/unadjusted_subplot_anova_stats[i,11]
    unadjusted_plot_anova_stats[i,12] <- unadjusted_plot_anova_stats[i,11]/unadjusted_subplot_anova_stats[i,11]
    

    
    #unadjusted_subplot_anova_stats[i,12] <- rxc_significance
    #unadjusted_plot_anova_stats[i,13] <- rxc_significance

    
  }
  
   
  anova_stats <-rbind(unadjusted_subplot_anova_stats,unadjusted_plot_anova_stats)
  anova_stats <- anova_stats[with(anova_stats , order(Year, Location, Type, Df)),]
  #anova_stats <- anova_stats[with(anova_stats , order(Type, Location, Year)),]
  
  GetMethod <- function(RS,CS,RXCS,RF,CF,RXFC) 
    {

    

    # determine method
    method='Method 1'
    if ((grepl( '*', RS, fixed = TRUE) && grepl( '*', CS, fixed = TRUE)) && grepl( 'ns', RXCS, fixed = TRUE)) 
      {
      method <-'Method 1'
      adj_method <-1
      } 
    else if ((grepl( 'ns', RS, fixed = TRUE)  && grepl( 'NS', CS, fixed = TRUE) ) && grepl( '*', RXCS, fixed = TRUE) ) 
     {
      method = 'Method 3'
      adj_method <-3
     } 
    else if ((grepl( 'ns', RS, fixed = TRUE)  && grepl( 'ns', CS, fixed = TRUE) ) && grepl( 'ns', RXCS, fixed = TRUE) ) 
     {
      method <-'Unnecessary'
      adj_method <-  0
      # According to ANOVA results, even row and column effects are not significant, still use Method 1
      #			method <-'Method 1'
      #			adj_method <-1
      } 
    else if ((grepl( '*', RS, fixed = TRUE)  && grepl( '*', CS, fixed = TRUE) ) && grepl( '*', RXCS, fixed = TRUE) ) 
      {
      method <-'Method 1'
      adj_method <- 1
      } 
    else if ((grepl( '*', RS, fixed = TRUE)  && grepl( 'ns', CS, fixed = TRUE) ) && grepl( '*', RXCS, fixed = TRUE) ) 
      {
      method <- 'Method 1 or Method 3 but'
      if (RXFC > CF && rxc_f > RF) 
        {
        method <- method + ' Method 3 is better'
        adj_method = 3
        } 
        else 
        {
        method <- method + ' Method 1 is better'
        adj_method = 1
        }
    }
    
    
    
    

    return(method)
  }
  
  traits <- unique(cp1$Trait)
  MethodsForTraits <- setNames(numeric(length(traits)), traits)
  
  for (i in 1:nrow(anova_stats)) #i=6
  {
    typeIs <- anova_stats[i,4]
    if (typeIs == "ERROR_WG") {
      df2Is<-anova_stats[i,5]
      rxc_fIs<-anova_stats[i,12]
      }
    if (typeIs == "ROW") {
      rowSig<-anova_stats[i,10]
      rowF<-anova_stats[i,8]
      }
    if (typeIs == "COL") {
      colSig<-anova_stats[i,10]
      colF<-anova_stats[i,8]
      }
    
    if (typeIs == "ERROR_RXC") {
      df1Is<-anova_stats[i,5]
      #rxc_fIs<-anova_stats[i,12]
      #anova_stats[i,3]
      fprob<- 1-pf(rxc_fIs, df1Is, df2Is,rowF,colF,rxc_fIs) #qf(rxc_fIs, df1Is, df2Is) 
      

       if (fprob <= 0.01) { 
         rxc_significance = "**"
       } else if (fprob <= 0.05) {
        rxc_significance = "*"
       } else {
        rxc_significance = "ns"
       }
      anova_stats[i,14]<-fprob
      methodVal<-GetMethod(rowSig,colSig,rxc_significance)
      anova_stats[i,15] <-methodVal
      MethodsForTraits[anova_stats[i,3]]<-methodVal
     }
    }
    

  

  
  write.csv(anova_stats , file = paste(DataFileDircsv,"MAD_ANOVA_result_summary",".csv"),row.names=FALSE)

###################adjusting_testlines_data_M1_M3
  
  #sub adjusting_testlines_data_M1_M3 
  p_reg_coes <- rcoeffs
  data_file <- datap1
  #adjusting_testlines_data_M1_M3
  #data_file
  #p_row_means2
  #p_col_means2
  #p_total_means
  #p_row_col_values
  #p_reg_coes
  data_file$adjmth1 <- 0  
  data_file$adjmth3 <- 0
   
  
  
  
  #####MTH 1 and MTH3 ###################################################################################3
  
   for(i in 1:nrow(data_file)) 
   {
     #i=9
     idis<-data_file[i,"IDX"]
     ids<-p_reg_coes["ListID"]
     idxpregcoes<-which(ids == idis)
     
     idisr<-data_file[i,"IDROW"]
     idr<-p_row_means2["IDROW"]
     idxrow<-which(idr == idisr)
     
     idisc<-data_file[i,"IDCOL"]
     idc<-p_col_means2["IDCOL"]
     idxcol<-which(idc == idisc)
     
     idisrc<-data_file[i,"IDROWCOL"]
     idrc<-p_row_col_values["IDROWCOL"]
     idxrowcol<-which(idrc == idisrc)
     
     if (data_file[i,12] == '.' || data_file[i,12] == 0) 
     {
      data_file[i,18] <- data_file[i,12]
      data_file[i,19] <- data_file[i,12]
     }
     else
     {
       data_file[i,18] = data_file[i,12] - p_row_means2[idxrow,"Value"] - p_col_means2[idxcol,"Value"] + 2*p_total_means[idxpregcoes,"Value"]
       xij <- p_row_col_values[idxrowcol,"Value"]
       if (!xij || is.na(xij)) 
       {
         xij <- (p_row_means2[idxrow,"Value"] + p_col_means2[idxcol,"Value"])/2;
       }

       
       
       
       if (p_reg_coes[idxpregcoes,"op"] && p_total_means[idxpregcoes,"Value"]) {
         data_file[i,19] <- data_file[i,12] - p_reg_coes[idxpregcoes,"op"] * (xij - p_total_means[idxpregcoes,"Value"]);
       }
       
       
       if (data_file[i,18] < 0 || is.na(data_file[i,18])) {data_file[i,18] <- 0}
       if (data_file[i,19] < 0 || is.na(data_file[i,19])) {data_file[i,19] <- 0}
     }
     
     
   }
  
  write.csv(data_file , file = paste(DataFileDircsv,"data_file_mth1_mth3",".csv"),row.names=FALSE)
  #####MTH 1 and MTH3 ###################################################################################3 
  
  
  
  
  
  
  
  
  #data_file$adjmth13x <- (data_file$adjmth3*100+0.5)/100
  #write.csv(data_file, file = paste(DataFileDircsv,"MAD_adjusted_all_data_Method1and3.csv",sep=""),row.names=FALSE)
  
  
  data_file$adjmth13 <- 0
  #####MTH1+3 ###################################################################################
  
  for(i in 1:nrow(data_file)) 
  {
    #i=9
    idis<-data_file[i,"IDX"]
    ids<-p_reg_coes["ListID"]
    idxpregcoes<-which(ids == idis)
    
    idisr<-data_file[i,"IDROW"]
    idr<-p_row_means2["IDROW"]
    idxrow<-which(idr == idisr)
    
    idisc<-data_file[i,"IDCOL"]
    idc<-p_col_means2["IDCOL"]
    idxcol<-which(idc == idisc)
    
    idisrc<-data_file[i,"IDROWCOL"]
    idrc<-p_row_col_values["IDROWCOL"]
    idxrowcol<-which(idrc == idisrc)
    
    if (data_file[i,18] == '.' || data_file[i,18] == 0) 
    {
      data_file[i,20] <- data_file[i,18]
    }
    else
    {
      if (p_reg_coes[idxpregcoes,"op"] && p_total_means[idxpregcoes,"Value"]) {
        data_file[i,20] <- data_file[i,18] - p_reg_coes[idxpregcoes,"op"] * (xij - p_total_means[idxpregcoes,"Value"]);
      }
      
      if (data_file[i,20] < 0 || is.na(data_file[i,20])) {data_file[i,20] <- 0}
    }
    
    
  }
  

  
  
  write.csv(data_file , file = paste(DataFileDircsv,"data_file_mth1+3",".csv"),row.names=FALSE)
  
  
  
  data_fileCRTL <- data_file
  data_fileCRTL$CTRLGENOTYPE <- 0 #create an environment attribute in the dataframe
  
  
  #if (data_fileCRTL$Cp == 1) 
  #{
  #  data_fileCRTL$CTRLGENOTYPE <- 1;
  #}
  #if (data_fileCRTL$Csp > 0) 
  #{
  #  data_fileCRTL$CTRLGENOTYPE <- 1;
  #}
  
  data_fileCRTL <- transform(data_fileCRTL, CTRLGENOTYPE= ifelse(Cp==1 | Csp > 0, 1, 0))
  
  write.csv(data_fileCRTL , file = paste(DataFileDircsv,"data_file_mth1+3+crtl",".csv"),row.names=FALSE)
  
  data_fileCRTLFilter<-subset(data_fileCRTL, CTRLGENOTYPE==1)
  
  newdata <- data_fileCRTLFilter[c(16,12,18:20)]
  
  dput(head(newdata))
  
  
  aggregate(Value~IDX,data=newdata,FUN=mean)
  meansq <- function(x) sum((x-mean(x))^2)/(length(x)-1)
  msdata <- newdata %>% group_by(IDX) %>% summarize(across(c(Value, adjmth1, adjmth3), meansq))
  write.csv(msdata , file = paste(DataFileDircsv,"ms",".csv"),row.names=FALSE)
  
  
  #####MTH1+3 ###################################################################################
  
  
  
  
  
  
  
  
  
  
  
  
  
  adjdata <- data_file
  ctldata <- adjdata
  
  ctldata <- filter(ctldata,(ctldata$Cp == 1) | (ctldata$Csp > 0))
  
  
  
#####################read_adjusted_data
  
  
  data_file_Adjustedm1m3 <- data_file
  data_file_Adjustedm1m3$control_genotypes <- 0
#  data_file_Adjustedm1m3$control_genotypes <- 0
  #datap1cg <- filter(data, (data$Csp > 0) | (data$Cp == 1))
  data_file_Adjustedm1m3FilteredCp <- filter(data_file_Adjustedm1m3,(data_file_Adjustedm1m3$Cp == 1))
  data_file_Adjustedm1m3FilteredCpTestPlots <- filter(data_file_Adjustedm1m3,(data_file_Adjustedm1m3$Cp != 1))
  #data_file_Adjustedm1m3FilteredCp <- filter(data_file_Adjustedm1m3,(data_file_Adjustedm1m3$Csp > 0) | (data_file_Adjustedm1m3$Cp == 1))
  #data_file_Adjustedm1m3FilteredCp <- filter(data_file_Adjustedm1m3,(data_file_Adjustedm1m3$Value != 0) )
  ctlgeno<-filter(data_file_Adjustedm1m3FilteredCp,(data_file_Adjustedm1m3FilteredCp$Cp == 1 | data_file_Adjustedm1m3FilteredCp$Csp > 0))
# write.csv(data_file_Adjustedm1m3FilteredCp, file = paste(txtdir,"data_file_Adjustedm1m3FilteredCp.csv",sep=""),row.names=FALSE)
  
  read_adjusted_data <- aggregate( adjmth1 ~ IDROW , data_file_Adjustedm1m3FilteredCp, mean )#row_means
  col_meansAM1 <- aggregate( adjmth1 ~ IDCOL , data_file_Adjustedm1m3FilteredCp, mean )#col_means
  total_meansAM1 <- aggregate( adjmth1 ~ IDX , data_file_Adjustedm1m3FilteredCp, mean )#total_means
  row_col_meansAM1 <- aggregate( adjmth1 ~ IDROWCOL , data_file_Adjustedm1m3FilteredCp, mean )#total_means
  row_col_testplots_meansAM1 <- aggregate( adjmth1 ~ IDROWCOL , data_file_Adjustedm1m3FilteredCpTestPlots, mean )#total_means
  row_meansAM3 <- aggregate( adjmth3 ~ IDROW , data_file_Adjustedm1m3FilteredCp, mean )#row_means
  col_meansAM3 <- aggregate( adjmth3 ~ IDCOL , data_file_Adjustedm1m3FilteredCp, mean )#col_means
  total_meansAM3 <- aggregate( adjmth3 ~ IDX , data_file_Adjustedm1m3FilteredCp, mean )#total_means
  row_col_meansAM3 <- aggregate( adjmth3 ~ IDROWCOL , data_file_Adjustedm1m3FilteredCp, mean )#total_means
  row_col_testplots_meansAM3 <- aggregate( adjmth3 ~ IDROWCOL , data_file_Adjustedm1m3FilteredCpTestPlots, mean )#total_means

  
  #Define the file name that will be deleted
  fn <- paste(DataFileDircsv,"read_adjusted_data",".xlsx",sep="")
  #Check its existence
  if (file.exists(fn)) {
    #Delete file if it exists
    file.remove(fn)
  }
  
  library(xlsx)
  write.xlsx(read_adjusted_data, file=fn, sheetName="read_adjusted_data", row.names=FALSE)
  write.xlsx(col_meansAM1, file=fn, sheetName="col_meansAM1", append=TRUE, row.names=FALSE)
  write.xlsx(row_col_meansAM1, file=fn, sheetName="row_col_meansAM1", append=TRUE, row.names=FALSE)
  write.xlsx(total_meansAM1, file=fn, sheetName="total_meansAM1", append=TRUE, row.names=FALSE)
  write.xlsx(row_col_testplots_meansAM1, file=fn, sheetName="row_col_testplots_meansAM1", append=TRUE, row.names=FALSE)
  write.xlsx(row_meansAM3, file=fn, sheetName="row_meansAM3", append=TRUE, row.names=FALSE)
  write.xlsx(col_meansAM3, file=fn, sheetName="col_meansAM3", append=TRUE, row.names=FALSE)
  write.xlsx(total_meansAM3, file=fn, sheetName="total_meansAM3", append=TRUE, row.names=FALSE)
  write.xlsx(row_col_meansAM3, file=fn, sheetName="row_col_meansAM3", append=TRUE, row.names=FALSE)
  write.xlsx(row_col_testplots_meansAM3, file=fn, sheetName="row_col_testplots_meansAM3", append=TRUE, row.names=FALSE)
     
  subplot_anova_file_data<-unadjusted_subplot_anova_stats
  
#ROBM
#####################################################################################################################################################

  subplot_ERROR <- read.csv(text="id,df,f,p", colClasses = c("integer", "integer", "numeric","numeric"))
  subplot_BLOCK <- read.csv(text="id,df,f,p", colClasses = c("integer", "integer", "numeric","numeric"))
  
  ctldata$IDXGT <- paste(ctldata$IDX,"_",ctldata$Genotype,sep="")
  control_values <- select(ctldata, c('IDXGT', 'IDX', 'Genotype', 'Value'))
  control_values_genotype <- unique(control_values$Genotype)
  control_values_m1 <- select(ctldata, c('IDXGT', 'IDX', 'Genotype', 'adjmth1'))
  control_values_genotype_m1 <- unique(control_values_m1$Genotype)
  control_values_m3 <- select(ctldata, c('IDXGT', 'IDX', 'Genotype', 'adjmth3'))
  control_values_genotype_m3 <- unique(control_values_m3$Genotype)
  
  
  control_values_m13 <- select(ctldata, c('IDXGT', 'IDX', 'Genotype', 'adjmth1'))
  #colnames(control_values_m13)[4] <- "adjmth13"
  control_values_m13$adjmth13  <- (control_values_m13$adjmth1*100+0.5)/100
  control_values_genotype_m13 <- unique(control_values_m3$Genotype)
  
  write.csv(ctldata, file = paste(DataFileDir,"ctldata",".csv",sep=""),row.names=FALSE)
  write.csv(ctldata, file = paste(DataFileDircsv,"ctldata",".csv",sep=""),row.names=FALSE)  
  
  vvv <- unique(ctldata$IDXGT)
  
  for(i in 1:nrow(subplot_anova_file_data)) #i=1
  {
    if (subplot_anova_file_data[i,4] == "ERROR_WG")
    {
      subplot_ERROR[i,1] <- paste(subplot_anova_file_data[i,1],subplot_anova_file_data[i,2],subplot_anova_file_data[i,3],sep="-")
      subplot_ERROR[i,2] <- subplot_anova_file_data[i,5]
      subplot_ERROR[i,3] <- subplot_anova_file_data[i,8]
      subplot_ERROR[i,4] <- subplot_anova_file_data[i,9]
    }
    else if ((subplot_anova_file_data[i,4] == "Whole_Plot"))
    {
      subplot_BLOCK[i,1] <- paste(subplot_anova_file_data[i,1],subplot_anova_file_data[i,2],subplot_anova_file_data[i,3],sep="-")
      subplot_BLOCK[i,2] <- subplot_anova_file_data[i,5]
      subplot_BLOCK[i,3] <- subplot_anova_file_data[i,8]
      subplot_BLOCK[i,4] <- subplot_anova_file_data[i,9]
    }
    else if ((subplot_anova_file_data[i,4] == "Genotype"))
    {
      subplot_BLOCK[i,1] <- paste(subplot_anova_file_data[i,1],subplot_anova_file_data[i,2],subplot_anova_file_data[i,3],sep="-")
      subplot_BLOCK[i,2] <- subplot_anova_file_data[i,5]
      subplot_BLOCK[i,3] <- subplot_anova_file_data[i,8]
      subplot_BLOCK[i,4] <- subplot_anova_file_data[i,9]
    }
    else
    {
      print  (paste("ERROR",subplot_anova_file_data[i,1],subplot_anova_file_data[i,2],subplot_anova_file_data[i,3],subplot_anova_file_data[i,4],sep="-"))
    }
  }
#####################################################################################################################################################    
  
  

  plot_anova_file_data<-unadjusted_plot_anova_stats 

  
  #debug
  write.csv(p_row_col_values2, file = paste(DataFileDir,"DEBUG_p_row_col_values2",".csv",sep=""),row.names=FALSE)
  write.csv(p_row_col_values , file = paste(DataFileDir,"DEBUG_p_row_col_values",".csv",sep=""),row.names=FALSE)
  #debug
  
  row_ids <- datap1$IDROW#row_ids
  col_ids <- datap1$IDCOL#col_ids  
  control_genotypes <- unique(datap1cg$Genotype) 
  
  

  
} else {
 file.remove(DataFilenameConverted)
  message("headers should be:Record	Plot	Row	Column	Cp	Csp	Entry	Year	Location	Genotype")
}


