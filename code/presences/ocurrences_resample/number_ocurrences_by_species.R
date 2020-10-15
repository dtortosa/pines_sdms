#Code that calculates the number of ocurrences by species and plot them as histogram

#####################################################################
########Calculation of the number of ocurrences per species##########
#####################################################################

###list ocurrences
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)  
#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false
list_species = cbind.data.frame(rep("Pinus", length(epithet_species_list)), epithet_species_list)
colnames(list_species) <- c("genus", "specific_epithet")
str(list_species)
head(list_species) #load data frame with genus and especific_epithet of each species. 

###loop for extract ocurrences number
number_ocurrences = NULL
for (i in list_species$specific_epithet){ #for each species
    ocurrences = read.csv(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(i, "final.presences.csv", sep="_"), sep="/")) #read the ocurrences
    number_ocurrences = append(number_ocurrences, nrow(ocurrences)) # calculate the number of rows and save it in empty vector number_ocurrences
}

#convert to a data frame
number_ocurrences = as.data.frame(number_ocurrences)

#bind name species.
number_ocurrences = cbind(list_species$specific_epithet, number_ocurrences)
names(number_ocurrences)[1] = "species"

#order from low to high number of ocurrences
number_ocurrences = number_ocurrences[with(number_ocurrences, order(number_ocurrences)),]

#test
nrow(number_ocurrences) == nrow(list_species) #one number for each species
str(number_ocurrences)
summary(number_ocurrences)
head(number_ocurrences)

#plot the result
plot(density((number_ocurrences$number_ocurrences)))

#extract the quantiles
quantile(number_ocurrences$number_ocurrences, c(0, 0.10, 0.20, 0.5, 0.9))

#write
write.csv(number_ocurrences, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/ocurrences_per_species.csv", row.names=FALSE)