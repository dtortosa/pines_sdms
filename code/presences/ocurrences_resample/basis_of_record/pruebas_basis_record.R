####Code for comprobate what is the number of fossil ocurrences

#########################################
##########Load species data##############
#########################################

########Tree species#########
#species by species using the bianca's species of the tree and gbif function of dismo
list_ocurrences = list.files(path="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species", pattern=".csv", full.names=TRUE)
length(list_ocurrences) #Load the path of each ocurrence file

#separated those with and without basis of record
record_table_1 = NULL
record_table_2 = NULL

for (i in 1:length(list_ocurrences)){
    table = read.csv(list_ocurrences[i], header=TRUE)
    if (length(table$basisOfRecord) > 0){ #If there is a column of basisOfRecord
        record_table_1 = rbind(record_table_1, table[,c("basisOfRecord", "scientificName", "lat", "lon")])
    } else { #if not
        record_table_2 = rbind(record_table_2, table[,c("scientificName", "lat", "lon")])
    }
}

str(record_table_1) 
str(record_table_2) #all species have a column of basisOfRecord

#drop rows without lon or without lat
record_table = record_table_1[!(is.na(record_table_1$lon) | is.na(record_table_1$lat)),]
str(record_table)
summary(record_table)

#check how many species with lat and long, have at the same time as value of basisOfRecord FOSSIL_SPECIMEN
nrow(record_table[record_table$basisOfRecord=="FOSSIL_SPECIMEN",]) #93 ocurrences 
record_table[record_table$basisOfRecord=="FOSSIL_SPECIMEN",]$scientificName #name of species from wich came these ocurrences. Most of fossiles from sylvestris and pseudo strobus. 

write.csv(record_table_1, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_resample/basis_of_record/basis_of_record_table.csv")

save.image("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/rdata/pruebas_basis_record.RData")