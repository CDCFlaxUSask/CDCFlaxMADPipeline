require(sasLM)

setwd("D:\\PlantSCI\\CDCFlaxMADPipeline")
d1 = read.csv(" MDA_plot_controls .csv")


d1 = af(d1, c("Year", "Location", "Trait")) 
f1 = Value ~ factor(Row) + factor(Column)
x4<-BY(GLM, f1, d1, By="Trait")

d2 = read.csv(" MDA_plot_subplot_controls .csv")
head(d2)


d2 = af(d2, c("Year", "Location", "Trait")) 
f2 = Value.x ~ factor(Whole_plot) + factor(Genotype.x)
x5<-BY(GLM, f2, d2, By="Trait")



setwd("D:\\PlantSCI\\CDCFlaxMADPipeline")
cp1 = read.csv("cp1s.csv")




#m1_row_col_values <-paste(cp1$year,"-",cp1$Location,"-",cp1$Trait),cp1$row,cp1$Column,cp1$Value#row_col_values


