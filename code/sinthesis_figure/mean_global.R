##########################################
######## CHANGE MEDIAN FOR MEAN ##########
##########################################

#load the data of suitability v2
differ_percent = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/differ_phylo_inside_v2.csv", sep=",", header=TRUE)
str(differ_percent)
head(differ_percent)

#calculate mean of all the columns except the first one with species name. Also remove the two first rows. These are the median and sd previously calculated
mean_values = cbind.data.frame("global mean", rbind.data.frame(apply(differ_percent[-c(1:2),-1], 2, mean)))
sd_values = cbind.data.frame("global standard deviation", rbind.data.frame(apply(differ_percent[-c(1:2),-1], 2, sd)))
names(mean_values) <- c("selected_species", colnames(differ_percent[,-1]))
names(sd_values) <- c("selected_species", colnames(differ_percent[,-1]))

#bind to the data.frame
differ_percent2 = rbind.data.frame(mean_values, sd_values, differ_percent)
str(differ_percent)

#check that the sd is the same in both dataframe. That row is calculated in the same way
round(differ_percent2[2,-1], 4) == round(differ_percent[2,-1], 4)

#check correlation between median and mean
differ_percent2[which(differ_percent2$selected_species == "global mean"),-1]
differ_percent[which(differ_percent2$selected_species == "global median"),-1]

#take a look
str(differ_percent2)
head(differ_percent2)

#save it
write.table(differ_percent2, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/differ_phylo_inside_v3.csv", col.names = TRUE, row.names = FALSE, sep=",")