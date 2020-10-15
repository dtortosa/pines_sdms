#code for making taxonomic comprobations about disagreements between taxonomy of Bianca and taxonomy of Gbif

#required package
require(dismo)

##############################################################
################ Taxonomic disagreements #####################
##############################################################
#####Taxonomic cleaning of the presences downlodaded from gbif of the entire genus. We have to solve taxonomic disagreements between taxonomy of Bianca and taxonomy of GBIF. 

###Load ocurrences directly downloades from gbif
pine_gbif = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/gbif_species/pines.csv", header=TRUE, sep="\t")
str(pine_gbif)
head(pine_gbif)
unique(pine_gbif$species) #154 species. There is 155 levels, but one of them is empty, there is no data from species neither speciesKey 

###Load list of pine species of the bianca's tree
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

#remove jaliscana and tecunumanii. These species were included in the new phylogeny of Bianca, but during my stay we don't have, neither phylogeny nor maps. These species were not included.
if("jaliscana" %in% epithet_species_list |  "tecunumanii" %in% epithet_species_list){
    epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("jaliscana"))]
    epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c( "tecunumanii"))]    
}
c("jaliscana", "tecunumanii") %in% epithet_species_list

#we add discolor because this species was used in the beginning
pos_discolor=28
epithet_species_list = c(epithet_species_list[1:(pos_discolor-1)], "discolor", epithet_species_list[pos_discolor:length(epithet_species_list)])
length(epithet_species_list) == 113
epithet_species_list[pos_discolor] == "discolor"

#load data frame with genus and especific_epithet of each species. 
list_species = cbind.data.frame(rep("Pinus", length(epithet_species_list)), epithet_species_list)
colnames(list_species) <- c("genus", "specific_epithet")
str(list_species)
head(list_species) 

#test with Pinus uncinata: 
unique(pine_gbif[pine_gbif$species== "Pinus mugo",]$species) #only mugo
unique(pine_gbif[pine_gbif$species== "Pinus mugo",]$specieskey) #there is two keys
nrow(pine_gbif[pine_gbif$species== "Pinus mugo",]) #only 7182 when gbif webpage indicate 21.882
nrow(pine_gbif[pine_gbif$specieskey== 5285385,]) #if we use the id of gbif webpage, there is more data, but don't reach 21.882
###conlcusion: DON'T use data directly download from the wab page. We will use the data download for my loop, species by species, EXECPT for pinus mugo, but I am sure that downloading the only the data of mugo from the webpage don't have fails. 

###Subset the species of gbif data set not included in the Bianca's tree and viceversa
species_out_bianca = unique(pine_gbif$species)[!unique(pine_gbif$species) %in% paste(list_species$genus, list_species$specific_epithet, sep=" ")] #Gbif's species that are considered as infraspecies or outside of Pinus according to Bianca's tree.  
species_out_gbif = paste(list_species$genus, list_species$specific_epithet, sep=" ")[!paste(list_species$genus, list_species$specific_epithet, sep=" ")  %in%  unique(pine_gbif$species)] #Bianca's tree species that are considered as infraspecies according to gbif. 
    
    #These species are the group for which we have problems. If a species is not included in our taxonomic concept but does in gbif, this could be a set of data that of Gbif bad assigned. The same goes for species included in our taxonomic concept, but not in gbif. In contrast, if a species is included in both gbif and my taxonomic concept, and it's not related to any of the disgreements, we can be sure that all the data included inside the name of that species is correct, i.e. no data belonging to other species included.

    #EXAMPLE: Pinus maritima Aiton is a non-accepted synonimum of P. halepensis, while Pinus maritima Lam. is a a non-accepted synonimum of P. pinaster. When we extracted data from gbif with species = "Pinus pinaster*", we have ocurrences of Pinus maritima Lam. included. Data coming from the request "Pinus halepensis*" has ocurrences of Pinus maritima Aiton included. In addition, if you search in the webpage of gbif for Pinus maritima Aiton, you see as the species name indicated is P. halepensis, while for Pinus maritima Lam. is P. pinaster. Therefore, GBIF takes into account its own taxonomic concept and it is not possible that a occurrence is included in P. halepensis* and P. pinaster*. Thus, if the taxonomic concept of gbif and Bianca's phylogeney matches, we can be sure that all entities included in that accepeted name are included in the search of gbif (Accepted name*).

    #In case a species included in our taxonomic concept but no in gbif was included in P. halepensis, we would detect it with this approach because this species does not match between both taxonomic concepts. In this case, we should extract these occurrence sand merge with the species of our taxonomic concept. 

    #If a species is in gbif but not in our taxonomic concept, then we should loook for other taxon in which our taxonomic concept includes that species. 

    #We compare the species names of our taxonimic concept with those found in request of ocurrences for the whole genus Pinus in one time yo gbif. This will include some species aceppted by gibf and not for us, but also some species that are not accepted for gbif even is included in the whole dataset of Pinus. this is the case of P. maritima. This species appears as a species in the search of the whole genus, but when we search species by species is included in halepensis or pinaster (for maritima Aiton and maritima Lam. respectively). We have no problem with these occurrences becasue we search species by species, we only use the list of species from the global pinus request to select potential disagreementes between taxonomic concepts. Each case is then evaluated separatelly using info from The plant list, Farjon's maps, gymnosperm database and gbif. 

    #Another example supporting that gbif search works. Pinus longifolia Salisb. is a non acepeted name of P. palustris, while Pinus longifolia Roxb. ex Lamb. is the non-accepted name of P. roxburgii. The specific search with the gbif function of the dismo package retreived data of Pinus longifolia Salisb. in palustris and data of Pinus longifolia Roxb. ex Lamb. This is correct.

    #Another one: Pinus tuberculata Gordon is included in P attenuata, while Pinus tuberculata D. Don is included in P. radiata. This is exactly correct. The same for Pinus muricata var. cedrosensis J. T. Howell, which is included in P. radiata, but not in P. muricata (See the plant list)

    #We could have problems if you different entities are accepted in Gbif, for example P. mugo and P. uncinata. However, in that case we would have one entity included in our taxonomic concept but not the other one, so we would detect it.

    #An aditional case: A species is considered as species for some authors but not for others. We don't have these species in the phylogeny and GBIF conisered as an subsespecies. This is the case of P. eremitana. This entity is considered as an species by Businsky, but Gbif included in P. fenzeliana (like other sources "http://www.catalogueoflife.org/col/details/species/id/a8bd6685755553a662215a0291cec17b"). If we don't have this entity (P. eremitana) as an distinct species and thus we don't have additional information from the taxonomic concept of Gallier paper, we follow the taxonomic concept og GBIf, which works very well. In addition, if any enetiy incuded by GBIF was very distant, remember we have our buffer that remove all ocurrences outside of it. This is the case of P. eremitana, its only ocurrence is removed, because it's outside the natural distribution.

    #An example of how GBif wors so well. Pinus maestrensis is included in cubensis for some authors, for others in occidentails and other considered as an species. GBIF said -> it belongs to cubensis. You can find occurrences of maestrensis in cubensis but not in occidentalis. It doest matter the taxonomic discussion, the decision of gibif is considered for all datasets. 

##############################################################
#############Cleaning presences from dismo loop###############
##############################################################
#####Taxonomic cleaning of the presences downlodaded from gbif function of dismo using a loop

###Gbif's species that are considered as infraspecies or outside of Pinus according to Bianca's tree.###### 
species_out_bianca 

###Pinus uncinata: 
#Nick include this species inside of pinus mugo, because there is not enough distance to consider them as different species, according to molecular markers. For further deatils see M. Heuertz, et al., Journal of Biogeography 37, 541 (2010). and I. Monteleone, D. Ferrazzini, P. Belletti, Silva Fennica 40, 391 (2006). Another interesting paper could be "Taxonomic revision of the Pinus mugo complex and P. × rhaetica (P. mugo × sylvestris) "
#If gbif considered as independnt species, its ocurrences are outside of Pinus mugo. In fact, there is not ocurrences of pinus mugo in SN according to bif, but yes of Pinus uncinata. 
pinus_uncinata = gbif(genus="Pinus", species="Pinus uncinata*", ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, ntries=5, nrecs=300, start=1, end=Inf) 
nrow(pinus_uncinata) #similar number to webpage (3.497)
unique(pinus_uncinata$speciesKey) #similar key to the webpage (5285332), a only one 
unique(pinus_uncinata$species) #only pinus uncinata
unique(pinus_uncinata$scientificName) #There is Pinus mugo subsp. uncinata, P. mugo var rostrata...

#load mugo data (directly downloaded from gbf webpage)
pinus_mugo = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_mugo.csv", header=TRUE)
nrow(pinus_mugo) #same number of ocurrences than number of ocurrences in the webpage (21893)
unique(pinus_mugo$specieskey) #a lot of keys 
unique(pinus_mugo$species) #only pinus mugo
unique(pinus_mugo$scientificname) #There is Pinus mugo subsp. uncinata (Ramond) Domin, Pinus uncinata var. rotundata (Link) Antoine, etc... It is to say, there is ocurrences of uncinata included here, but not all, because for example there are not ocurrences in SN or central Spain (http://www.gbif.org/species/5285385), on the contrary of P.uncinata (http://www.gbif.org/species/5285332).

#This is not usual in GBIf because, in all the checks I have done, occurrences of the accepted species in GBIF are downlodaded only. Maybe could be cause because I donwlodaded these occurrences from the webpage instead of the function "gbif" of dismo.

#Solution: Merge both dataset for include ALL ocurrences of P.uncinata. There is no problem in relation to duplicated ocurrences, because we are going to delete duplicated, and make the resampling of ocurrences. In addition, points extranges in water are eliminaited because all points without soil data are removed. 

#see both structures
str(pinus_mugo)
str(pinus_uncinata)

#make the merge
pinus_mugo_clean = pinus_mugo[,c("basisofrecord", "scientificname", "specieskey",  "decimallongitude", "decimallatitude", "year", "coordinateprecision", "coordinateuncertaintyinmeters")]
pinus_uncinata_clean = pinus_uncinata[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinatePrecision", "coordinateUncertaintyInMeters")] #select interesting columns
names(pinus_mugo_clean) <- c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinatePrecision", "coordinateUncertaintyInMeters")
names(pinus_uncinata_clean) <- c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinatePrecision", "coordinateUncertaintyInMeters") #change the names

#see structure
str(pinus_mugo_clean)
str(pinus_uncinata_clean)

#make rbind. 
pinus_mugo_final = rbind(pinus_mugo_clean, pinus_uncinata_clean)
str(pinus_mugo_final)
nrow(pinus_mugo_final) == nrow(pinus_mugo_clean) + nrow(pinus_uncinata_clean) #same number of rows
ncol(pinus_mugo_final) == ncol(pinus_mugo_clean)
ncol(pinus_mugo_final) == ncol(pinus_uncinata_clean) #same number of columns

#write the csv ONLY ONE TIME 
write.csv(pinus_mugo_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_mugo.csv", row.names = FALSE)

#plot points: Uncinata in SN. 
pinus_mugo_final = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_mugo.csv", sep=",", header=T)
bio1 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio1.asc")
plot(crop(bio1, c(-10,40,30,70)), col="gray")
points(pinus_mugo_final$lon, pinus_mugo_final$lat, cex=0.5)

#Pinus densithunbergii: Pinus × densithunbergii, first described by Homiki Uyeki in 1926 is the name given to the grex of Pinus densiflora and Pinus thunbergii. Common names include アイグロマツ (Aiguro-matsu) and アカクロマツ (Akakuro-matsu) in the Japanese language. Some publications will use the synonym, Pinus × densi-thunbergii (see http://conifersociety.org/conifers/conifer/pinus/x-densithunbergii/). It is a hybrid of two species, because of this we will not use it. Moreover, in gbif there is only 30 ocurrences. 

#Pinus rhaetica. It is a hybrid between Pinus sylvestris and Pinus mugo subps. uncinata (see Taxonomic revision of the Pinus mugo complex and P. x rhaetica (P. mugo x sylvestris) (Pinaceae) Christensen 1996). Only 44 ocurrences in gbif. We don't use it. 
pinus_sylvestris = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_sylvestris.csv", sep=",", header=T)
unique(pinus_mugo_final$scientificName)
unique(pinus_sylvestris$scientificName) #Pinus rhaetica not included in sylvestris or mugo.

#Pinus hugeli. There is no information about this species, it could be an error, and really would means Pinus hugelii, a synonymum of Pinus teocote according to gbif (http://www.gbif.org/species/5285233/synonyms). There is only two ocurrences in gbif. OUT. 
pinus_teocote = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_teocote.csv", sep=",", header=T)
unique(pinus_teocote$scientificName) #not included

#Pinus grossana. It is a fossil pine species (http://www.gbif.org/occurrence/search?taxon_key=5822776). OUT. 

#Pinus saturni. It is a fossil pine species (http://www.gbif.org/occurrence/search?taxon_key=5822777). OUT. 

#Pinus larix. It is a fossil pine species (http://www.gbif.org/occurrence/search?taxon_key=7458308) or a synonimun of Larix decidua? (http://www.theplantlist.org/tpl1.1/record/kew-2561645) OUT. 

#Pinus tecunumanii. It had been considered as a variety of pinus patula or pinus oocarpa but this has been refused by RAPD markers (furman et al 1997, analysis of genetic relatioships...). It has 348 ocurrences in gbif. It is diffrent species of the 113 of Bianca, thus we don't have buffer nor phylo tree. OUT
pinus_oocarpa = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_oocarpa.csv", header=TRUE)
pinus_patula = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_patula.csv", header=TRUE)
unique(pinus_oocarpa$scientificName)
unique(pinus_patula$scientificName) #This species is no included in oocpara or patual for gbif data. 

#pinus jaliscana. It is very close to Pinus oocarpa, Pinus herrerae and Pinus patula, but it seems that they are different species (Dvorak, W.S., A.P. Jordon, G.P. Hodge and J.L. Romero. 2000. Assessing evolutionary relationships of pines in the Oocarpae and Australes subsections using RAPD markers. New Forests 20(2):163-192.). Only 56 records. It is diffrent species of the 113 of Bianca, thus we don't have buffer. OUT 
pinus_oocarpa = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_oocarpa.csv", header=TRUE)
pinus_patula = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_patula.csv", header=TRUE)
pinus_herrerae = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_herrerae.csv", header=TRUE)
unique(pinus_oocarpa$scientificName)
unique(pinus_patula$scientificName) 
unique(pinus_herrerae$scientificName) #This species is no included in oocarpa or patula for gbif data. 

##Pinus luzmariae: This taxon was until recently known as Pinus oocarpa var. trifoliata (Farjon, A., and B.T. Styles. 1997. Pinus (Pinaceae). Flora Neotropica Monograph 75. New York, NY: The New York Botanical Garden.), but Pérez de la Rosa (1998) (Pérez de la Rosa, J.A. 1998. Promoción de una variedad de pino serotino mexicano a nivel de especie. Bol. Inst. Bot. Univ. Guadalajara 5: 127-135.) has made a case for its recognition as a species, finding additional characters besides the numbers of needles in a fascicle to distinguish it. The ecology is similar to oocarpa ("An Atlas of the World's Conifers: An Analysis of their Distribution, Biogeography, Diversity and Conservation Status"). The scattered distribution across a large part of Mexico and beyond, always within the very wide range of P. oocarpa, and the rarity of the trees with these characters in any location, seem to go against this, but it can be given the benefit of the doubt in the absence of evidence, e.g. based on DNA data, to the contrary. I can`t see evidence for separate, and there is only 42 ocurrences. 
#According to Farjon atlas, most points of this species are included inside our buffer of oocarpa
pinus_oocarpa = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_oocarpa.csv", header=TRUE)
unique(pinus_oocarpa$scientificName) #This species is no included in oocpara or patual for gbif data. 

#download
pinus_luzmariae = gbif(genus="Pinus", species="Pinus luzmariae*", ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, ntries=5, nrecs=300, start=1, end=Inf)
unique(pinus_luzmariae$scientificName) 
unique(pinus_luzmariae$speciesKey) #There is "Pinus oocarpa var. trifoliata Martinez" and "Pinus luzmariae Pérez de la Rosa". Both enttites fall inside luzmariae ("An Atlas of the World's Conifers: An Analysis of their Distribution, Biogeography, Diversity and Conservation Status")

#We are going to merge both dataset, because I don't have enough evidence to consider it as a different species.
pinus_oocarpa_clean = pinus_oocarpa[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_luzmariae_clean = pinus_luzmariae[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year")]
names(pinus_oocarpa_clean) <- c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")
names(pinus_luzmariae_clean) <- c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year")
pinus_oocarpa_clean$coordinatePrecision = rep(NA, times=nrow(pinus_oocarpa_clean))
pinus_luzmariae_clean$coordinatePrecision = rep(NA, times=nrow(pinus_luzmariae_clean))
pinus_luzmariae_clean$coordinateUncertaintyInMeters = rep(NA, times=nrow(pinus_luzmariae_clean))
pinus_oocarpa_final = rbind(pinus_oocarpa_clean, pinus_luzmariae_clean)
str(pinus_oocarpa_final)
write.csv(pinus_oocarpa_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_oocarpa.csv", row.names = FALSE)

#Pinus longifolia. Pinus longifolia Salisb is a synonimus of P. palustris, while Pinus longifolia Roxb. ex Lamb. is a synonimus of pinus_roxburghii.
#download longifolia
pinus_longifolia = gbif(genus="Pinus", species="Pinus longifolia*", ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, ntries=5, nrecs=300, start=1, end=Inf)
unique(pinus_longifolia$scientificName) 
unique(pinus_longifolia$speciesKey)
unique(pinus_longifolia$lon)
unique(pinus_longifolia$lat) #there is not ocurrences with coordiantes

#take a look inside pinus_roxburghii and pinus_palustris
pinus_roxburghii = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_roxburghii.csv", header=TRUE)
unique(pinus_roxburghii$scientificName) #Pinus longifolia Roxb. ex Lamb. is already included in pinus_roxburghii. 
pinus_palustris = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_palustris.csv", header=TRUE)
unique(pinus_palustris$scientificName) #Pinus longifolia Salisb is already included in this species. 

##Pinus attenuradiata: It is a hybrid between attenuata and radiata (http://plants.usda.gov/core/profile?symbol=piat2). OUT.
pinus_attenuata = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_attenuata.csv", header=TRUE)
pinus_radiata = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_radiata.csv", header=TRUE)
unique(pinus_attenuata$scientificName)
unique(pinus_radiata$scientificName) #not included in these species. OJO, Pinus muricata var. cedrosensis se incluye en radiata "https://www.fs.fed.us/database/feis/plants/tree/pinrad/all.html"
pinus_muricata = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_muricata.csv", header=TRUE)
unique(pinus_muricata$scientificName) #Pinus muricata D. Don var. cedrosensis not included

##Pinus georginae: It is a species recent discover by the la Rosa (Pinus georginae (Pinaceae), a new species from western Jalisco, Mexico, De la Rosa 2009).Only 16 ocurrences. Before it has bee included in P. praetermissa and genetically close to P. luzumariae.It is diffrent species of the 113 of Bianca, thus we don't have buffer. OUT
pinus_praetermissa = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_praetermissa.csv", header=TRUE)
unique(pinus_praetermissa$scientificName) #Pinus georginae not included
pinus_oocarpa_final = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_oocarpa.csv", header=TRUE)
unique(pinus_oocarpa_final$scientificName) #Pinus georginae not included in oocarpa, a species with data from luzumariae 


#Pinus abies. It seems that this is an error, there is Picea abies, but not pinus abies. Only 3 ocurrences. OUT. 

#Pinus neilreichiana. It is an hybrid (http://conifersociety.org/conifers/conifer/pinus/x-neilreichiana/). The natural hybrid between P. nigra subsp. nigra and Pinus sylvestris. No synonyms. Type: Austria, Niederösterreich, Baden, Grossau, along footpath to Pottenstein, Reichardt s.n., lectotype C (https://www.conifers.org/pi/Pinus_nigra.php) Only 8 ocurrences. OUT. 

#Pinus maritima. P. maritima Aiton and Pinus maritima Mill. belongs to halepenis, while Pinus maritima Lam. belongs to P. pinaster
pinus_halepensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_halepensis.csv", header=TRUE)
unique(pinus_halepensis$scientificName) #Pinus maritima. P. maritima Aiton and Pinus maritima Mill. included
pinus_pinaster = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_pinaster.csv", header=TRUE)
unique(pinus_pinaster$scientificName) #Pinus maritima Lam. included

#Pinus austriaca. It is a other name of Pinus nigra, not accepted. Only 2 georreferenciated points and we have 33560 points of Pinus nigra. OUT. 
pinus_nigra = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_nigra.csv", header=TRUE)
unique(pinus_nigra$scientificName)#Pinus austriaca Höss, Pinus nigra var. austriaca (Höss) Badoux, Pinus laricio subsp. austriaca (Höss) Endl. included

#Pinus montana. Probably a synonimun of Pinus cembra, sylvestris or mugo. Only 1 point with coordinates, and we have a lot of points for these species. OUT. 
pinus_cembra = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_cembra.csv", header=TRUE)
unique(pinus_cembra$scientificName)# Pinus montana Lam. included      

#Pinus hakkodensis. P. × hakkodensis Makino 1931 is the natural hybrid between P. parviflora var. pentaphylla and P. pumila. Synonym: P. pentaphylla Mayr var. hakkodensis (Makino) Kusaka 1954 (gymnosperm database). Only 8 points with coordinates. OUT.
pinus_parviflora = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_parviflora.csv", header=TRUE)
unique(pinus_parviflora$scientificName) #not included
pinus_pumila = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_pumila.csv", header=TRUE)
unique(pinus_pumila$scientificName) #not included


#Pinus henryi. This species has been treated as a variety of Pinus tabuliformis in the Flora of China, Vol. 4 (1999) and is indeed similar in many characters to that species. Its leaves are not wider than 1 mm and its seed cones are smaller, with scales that have only slightly raised apophyses. It was also classified as a variety of P. massoniana, from which it differs in its much shorter leaves, its uninodal growth of shoots, smaller, nearly globose-ovoid seed cones, and relatively shorter seed wings. I am not sure where included these ocurrences and they are only 30 points with lon/lat. Moreover, microsatellyte data prompt to genetical differentiation between theses species (High genetic differentiation in natural populations of Pinus henryi and Pinus tabuliformis as revealed by nuclear microsatellites, Zhan-Lin Liu et al-. 2012). I have considered as independent species, and thus OUT (no buffer).   
pinus_tabuliformis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_tabuliformis.csv", header=TRUE)
pinus_massoniana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_massoniana.csv", header=TRUE)
unique(pinus_tabuliformis$scientificName)
unique(pinus_massoniana$scientificName) #not included in these species

#Pinus khasya. A synonimun of Pinus kesiya Royle ex Gordon. Only 5 georeferenced points and we have 369 points of Pinus kesiya. OUT. 
pinus_kesiya = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_kesiya.csv", header=TRUE)
unique(pinus_kesiya$scientificName)# Pinus khasya Hook. f.

#Pinus schwerinii. Hybrid of Pinus strobus and P. Wallichiana. No georefenced points. OUT. 
pinus_strobus = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_strobus.csv", header=TRUE)
pinus_wallichiana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_wallichiana.csv", header=TRUE)
unique(pinus_strobus$scientificName)
unique(pinus_wallichiana$scientificName) #not included in these species 

#Pinus excelsa. Synonimun of P. Wallichiana. Only 6 points with coordinates of which only 2 are inside of the natural distribution area. Moreover we have 541 ocurrences of P. Wallichiana. OUT. 
    #OJO: En the plant list pone que es sinonimo de picea abies!!! http://www.theplantlist.org/tpl1.1/record/kew-2562069
pinus_wallichiana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_wallichiana.csv", header=TRUE)
unique(pinus_wallichiana$scientificName)# Pinus excelsa Wall. ex D. Don

#Pinus sinensis. Synonimun of P. massoniana. There is not points with coordinates. OUT. 
pinus_massoniana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_massoniana.csv", header=TRUE)
unique(pinus_massoniana$scientificName)# Pinus sinensis D. Don included

#Pinus insignis. Synonimun of P. radiata. Only 1 point with coordiantes. OUT. 
pinus_radiata = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_radiata.csv", header=TRUE)
unique(pinus_radiata$scientificName)# Pinus insignis Douglas ex Loudon and Pinus insignis var. macrocarpa Hartw. included

#Pinus pallasiana. Variety of P. nigra and any point woth coordinates. OUT. 
pinus_nigra = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_nigra.csv", header=TRUE)
unique(pinus_nigra$scientificName)#Pinus nigra var. pallasiana (Lamb.) C. K. Schneid., Pinus nigra var. pallasiana (Lamb.) Asch. & Graebn., Pinus pallasiana Lamb., Pinus laricio var. pallasiana (Lamb.) Endl., Pinus nigra subsp. pallasiana (Lamb.) Holmboe

#Pinus taihangshanensis. Synonim of P. tabuliformis ("http://www.theplantlist.org/tpl1.1/record/kew-2562319"). No coordinates data. OUT.
pinus_tabuliformis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_tabuliformis.csv", header=TRUE)
unique(pinus_tabuliformis$scientificName)#Not included

#Pinus holfordiana. Hybrid of Pinus ayacahuite and Pinus wallichiana (http://www.botanic.cam.ac.uk/Botanic/Plant.aspx?p=27&ix=37&pid=0&prcid=0&ppid=0). 2 coordinates data. OUT. 
pinus_ayacahuite = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_ayacahuite.csv", header=TRUE)
pinus_wallichiana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_wallichiana.csv", header=TRUE)
unique(pinus_ayacahuite$scientificName)
unique(pinus_wallichiana$scientificName) #not included in these species

#Pinus sondereggeri. Hybrid of Pinus taeda and P. palustris (https://en.wikipedia.org/wiki/Pinus_%C3%97_sondereggeri). Without georeferenced data. OUT. 
pinus_taeda = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_taeda.csv", header=TRUE)
pinus_palustris = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_palustris.csv", header=TRUE)
unique(pinus_taeda$scientificName)
unique(pinus_palustris$scientificName) #not included in these species

#Pinus longipedunculata. Variety of Pinus patula var. longipedunculata Loock ex Martínez. Without georeferenced data. OUT. 
pinus_patula = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_patula.csv", header=TRUE)
unique(pinus_patula$scientificName)#Pinus patula var. longipedunculata Loock ex Martínez and Pinus patula var. longipedunculata Loock included

#Pinus pithyusa. Variety of Pinus brutia var. pityusa (Steven) Silba. Without georeferenced data. OUT. 
pinus_brutia = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_brutia.csv", header=TRUE)
unique(pinus_brutia$scientificName)#Pinus pityusa Steven and Pinus pithyusa Fox-Strangw. included

#Pinus hwangshanensis is the Chinese representative of a group of three closely related taxa: Pinus taiwanensis of Taiwan, Pinus luchuensis of Japan, and Pinus hwangshanensis of mainland China. All three taxa are morphologically similar, but distinct, and they are here treated as separate species, although they could also be called subspecies of P. luchuensis. Only 4 georeferenced data. I don't know at which species belongs these ocurrences. OUT. 
pinus_taiwanensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_taiwanensis.csv", header=TRUE)
pinus_luchuensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_luchuensis.csv", header=TRUE)
unique(pinus_taiwanensis$scientificName)
unique(pinus_luchuensis$scientificName) #hwangshanensis is not included in this species

#Pinus rudis. Synonimun of Pinus hartwegii. Only 3 points with coordinates, and we have 1188 points of P.hartwegii. OUT.
pinus_hartwegii = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_hartwegii.csv", header=TRUE)
unique(pinus_hartwegii$scientificName)# Pinus rudis Endl. included

#Pinus australis. Synonimun of Pinus palustris. Without georeferenced data. OUT.
pinus_palustris = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_palustris.csv", header=TRUE)
unique(pinus_palustris$scientificName)# Pinus australis F. Michx. included

#Pinus inops. Synonimun of P. contorta or virginiana. Without georeferenced data. OUT.
pinus_clausa = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_contorta.csv", header=TRUE)
unique(pinus_clausa$scientificName)# Pinus inops Bong. included
pinus_virgniana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_virginiana.csv", header=TRUE)
unique(pinus_virgniana$scientificName)# Pinus inops Aiton included

#Pinus picea. Synonimum of Picea abies. Without georeferenced data. OUT.

#Pinus scotica. Variety of sylvestris. Without georeferenced data. OUT.
pinus_sylvestris = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_sylvestris.csv", header=TRUE)
unique(pinus_sylvestris$scientificName)# Pinus sylvestris var. scotica Beissn. included

#Pinus iddingsi. Fossil species. OUT. 

#Pinus peregrinus. Fossil species. OUT. 

#Pinus estevezi. Synonimun of Pinus pseudostrobus Lindl. Only 1 point with coordinates and we have 1330 point of P. pseudostrobus. OUT. 
pinus_pseudostrobus = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_pseudostrobus.csv", header=TRUE)
unique(pinus_pseudostrobus$scientificName)# Pinus estevezii (Martínez) J. P. Perry included

#Pinus tecumumani. Wrong of tecunumanii. No points with coordinates. OUT.

#Pinus menziesii. Wrong, is Picea sitchensis (Bong.) Carrière. OUT. 

#Pinus pentaphylla. A variety of P. parviflora. Only 2 points with coordinates, and only one in Japan (endemic of this country). We have 600 points of P. parviflora. OUT. 
pinus_parviflora = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_parviflora.csv", header=TRUE)
unique(pinus_parviflora$scientificName)# Pinus pentaphylla Mayr and Pinus parviflora var. pentaphylla (Mayr) A. Henryincluded

#Pinus parvifolia. It seems an error, the correct is P. parviflora. Only 5 points with coordinates and any one inside natural distribution (Japan), and parviflora tiene 600 datos. 
pinus_parviflora = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_parviflora.csv", header=TRUE)
nrow(pinus_parviflora)

#Pinus maximowiczii. It is from Picea.  Without georeferenced data. OUT.

#Pinus wangii. It is a endemism of China. 23 points with coordinates. It is diffrent species of the 113 of Bianca, thus we don't have buffer. OUT
unique(pinus_parviflora$scientificName)#not included

#Pinus divaricata. Synonimun of P. banksiana. Only 2 points with coordinates and all outisde of the natural distribution. Moreover, we have 1848 points of P. banksiana. OUT. 
pinus_banksiana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_banksiana.csv", header=TRUE)
unique(pinus_banksiana$scientificName)#Pinus divaricata (Aiton) Dum.-Cours. and Pinus divaricata (Aiton) Sudw. included

#Pinus pinsapo. Error, it is Abies pinsapo. Without georeferenced data. OUT.

#Pinus pygmaea. Pinus cembra var. pygmaea Loudon is synonimum of P. pumila. No points. OUT. 
pinus_pumila = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_pumila.csv", header=TRUE)
unique(pinus_pumila$scientificName)#Not included here

#Pinus finlaysoniana. Synonimum of Pinus merkusii Jungh. & de Vriese. No points with coordinates. OUT. 
pinus_merkusii = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_merkusii.csv", header=TRUE)
unique(pinus_merkusii$scientificName)#Not included here


###Bianca's tree species which are not considered as species in gbif###### 
species_out_gbif
#All Bianca's trees speceis have been download with the loop withput problems. However, some species are not downlodaded if we use the request data for the entire genus. This is due to whe you ask for Pinus in gbif only download the accepted speceis according to gbif. But, inside of these speceis are included species don't accepted by gbif but yes for Nick. For example: Pinus cooperi is not downloaded with Pinus request, but inside P. arizonica there ocurrences with scientificName of Pinus arizonica subsp. cooperi. It is to say, gbif consider Pinus cooperi as a subspecies of Pinus arizonica. However, Bianca`s tree has a P.cooperi, thus I included in the loop. When you take a look inside cooperi data.frame there is ocurrences with scientificName=P.cooperi but with the speciesCode of P. arizonica (5285094). Therefore, when we request for "Genus epithet_specific*", we obtain al ocurrences with the name of the species in scientificName although this species are not accepted by gbif and it is inside an accepted species. Conclusion: In Pinus_cooperi.csv are included all ocurrences of P.cooperi, BUT theses ocurrences are also included in Pinus_arizonica.csv. Therefore, we have to drop the presences of the non-gbif accepted species from the data.frame of the accepted species. 

#I checked which species move from a entity to another using the plant list, for exmaple: Pinus lutea C. E. Blanco ex Martínez is inside arizonica, but belongs to cooperii, even cooper is not in the name. We know this from The Plant List

#Pinus amamiana. This taxa is considered as speceis according to gbif (http://www.gbif.org/species/5285748). There is 6 ocurrences. No problem. 
pinus_amamiana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_amamiana.csv", header=TRUE)
unique(pinus_amamiana$speciesKey) 
unique(pinus_amamiana$scientificName)
    #[1] Pinus amamiana Koidz. -> amamiana                
    #[2] Pinus armandii var. amamiana (Koidz.) Hatus. -> amamiana

#Pinus chiapensis. It is included in Pinus strobus according to gbif. Section Cembra, the 5-needled white pines. Syn: P. strobus L. var. chiapensis Martínez 1940; P. strobus L. subsp. chiapensis (Martínez) E. Murray 1982. Its close relationship to P. strobus L., of Canada and the eastern U.S., has been substantiated using molecular as well as morphological and chemical evidence. Many authors continue to treat it as a subspecies or variety of P. strobus. The small morphological differences between the two taxa, discussed by Farjon and Styles (1997), can all be accounted for by habitat differences (the taxa are separated by a distance of no less than 2000 km, and by the presence/absence of seasonal frost). It seems likely that we are looking at a classic example of vicariant speciation, in which continued geographic separation between the two taxa would allow them to develop into clearly distinct species, while elimination of the separation (perhaps via a return to the conditions that existed during the last glacial maximum) might allow interbreeding and the disappearance of P. chiapensis as a distinct taxon. Molecular data have also established that these two taxa are closely related to the other Mexican five-needle white pines, P. ayacahuite, P. strobiformis and P. flexilis (Liston et al. 1999, Gernandt et al. 2005). They are DIFFERENTS. Source: "The Gymnosperm Database"
pinus_strobus = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_strobus.csv", header=TRUE)
unique(pinus_strobus$speciesKey) 
unique(pinus_strobus$scientificName) #it is inlcuded chiapensis here. 
    #[1] Pinus strobus L. -> strobus                                    
    #[2] Pinus strobus var. chiapensis Martinez -> chiapensis             
    #[3] Pinus strobus subsp. cumberlandensis Silba -> strobus          
    #[4] Pinus strobus var. strobus -> strobus                         
    #[5] Pinus chiapensis (Martínez) Andresen -> chiapensis               
    #[6] Pinus strobus subsp. chiapensis (Martínez) A.E.Murray -> chiapensis
    #[7] Pinus strobus subsp. chiapensis (Martínez) E. Murray -> chiapensis

pinus_chiapensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_chiapensis.csv", header=TRUE)
unique(pinus_chiapensis$speciesKey) #the same key of Pinus strobus 
unique(pinus_chiapensis$scientificName) #If you take a look on the distributon maps of ocurrences of P.strobus (http://www.gbif.org/species/5284982) and P. chiapensis (http://www.gbif.org/species/5284987). There are some ocurrences of P. strobus chiapensis  dataset don't included in P.chiapensis data set. Therefore, we have to drop these ocurrences from strobus and included in chiapensis. No problem with duplicated (the will be drop with te code later). 
    #Pinus chiapensis (Martínez) Andresen -> chiapensis

#drop ocurrences of chiapensis from strobus
pinus_strobus_clean = pinus_strobus[!(pinus_strobus$scientificName=="Pinus strobus var. chiapensis Martinez" | pinus_strobus$scientificName=="Pinus chiapensis (Martínez) Andresen" | pinus_strobus$scientificName=="Pinus strobus subsp. chiapensis (Martínez) A.E.Murray" | pinus_strobus$scientificName=="Pinus strobus subsp. chiapensis (Martínez) E. Murray"),] #subset by scientificName with chiapensis
unique(pinus_strobus_clean$scientificName)

#write it
write.csv(pinus_strobus_clean, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_strobus.csv", row.names = FALSE)

#include these ocurrences in chiapensis
pinus_chiapensis_strobus = pinus_strobus[(pinus_strobus$scientificName=="Pinus strobus var. chiapensis Martinez" | pinus_strobus$scientificName=="Pinus chiapensis (Martínez) Andresen" | pinus_strobus$scientificName=="Pinus strobus subsp. chiapensis (Martínez) A.E.Murray" | pinus_strobus$scientificName=="Pinus strobus subsp. chiapensis (Martínez) E. Murray"),]

#make the merge
pinus_chiapensis_clean = pinus_chiapensis[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year")]
pinus_chiapensis_strobus_clean = pinus_chiapensis_strobus[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinatePrecision", "coordinateUncertaintyInMeters")]
pinus_chiapensis_clean$coordinatePrecision = rep(NA, times=nrow(pinus_chiapensis_clean))
pinus_chiapensis_clean$coordinateUncertaintyInMeters = rep(NA, times=nrow(pinus_chiapensis_clean))
pinus_chiapensis_final = rbind(pinus_chiapensis_clean, pinus_chiapensis_strobus_clean)
str(pinus_chiapensis_final)
nrow(pinus_chiapensis_final) == nrow(pinus_chiapensis_clean) + nrow(pinus_chiapensis_strobus_clean)
unique(pinus_chiapensis_final$scientificName)

#write it
write.csv(pinus_chiapensis_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_chiapensis.csv", row.names = FALSE)

#Pinus cooperi. It is included in arizonica by gbif. Farjon and Styles (1997) assert that, as with the other varieties of P. arizonica, characters distinguishing the varieties vary continuously, thus treatment at the species rank is inappropriate. However, quite a few Mexican pine species are not clearly distinguishable from their close relatives, due to processes such as hybridization and ongoing species differentiation, so this argument is rejected by some Mexican botanists and the name Pinus cooperi is in fairly common use in the literature. It may also be appropriate to treat this taxon at the subspecies rank, and it has been so described, but the name is not widely used. . It is not clear, in the bibliography both terms are used, but there is not genetic information. I will use it as different species for using their buffer. Source: "The Gymnosperm Database" 
pinus_arizonica = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_arizonica.csv", header=TRUE)
unique(pinus_arizonica$speciesKey) 
unique(pinus_arizonica$scientificName) #it is inlcuded cooperi here. 
    #Entities included in pinus_arizonica
        #[1] Pinus arizonica Engelm. -> P. arizonica                           
        #[2] Pinus arizonica var. stormiae Martinez -> Pinus arizonica 
        #[3] Pinus arizonica var. arizonica -> P. arizonica                    
        #[4] Pinus arizonica var. cooperi (C. E. Blanco) Farjon -> cooperi 
        #[5] Pinus ponderosa var. arizonica (Engelm.) Shaw -> arizonica    
        #[6] Pinus cooperi var. ornelasii (Martínez) C. E. Blanco -> cooperi
        #[7] Pinus cooperi C. E. Blanco -> cooperi                      
        #[8] Pinus lutea C. E. Blanco ex Martínez -> cooperi             
        #[9] Pinus ponderosa var. stormiae (Martínez) Silba -> arizonica     
        #[10] Pinus ponderosa subsp. arizonica (Engelm.) E. Murray -> arizonica

pinus_cooperi = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_cooperi.csv", header=TRUE)
unique(pinus_cooperi$speciesKey) 
unique(pinus_cooperi$scientificName) #only cooperi
    #Pinus cooperi C. E. Blanco -> cooperi

#drop ocurrences of cooperi from arizonica
pinus_arizonica_clean = pinus_arizonica[!(pinus_arizonica$scientificName=="Pinus arizonica var. cooperi (C. E. Blanco) Farjon" | pinus_arizonica$scientificName=="Pinus cooperi var. ornelasii (Martínez) C. E. Blanco" | pinus_arizonica$scientificName=="Pinus cooperi C. E. Blanco" | pinus_arizonica$scientificName=="Pinus lutea C. E. Blanco ex Martínez"),]
unique(pinus_arizonica_clean$scientificName)

#write it
write.csv(pinus_arizonica_clean, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_arizonica.csv", row.names = FALSE)

#include these ocurrences in cooperi 
pinus_cooperi_arizonica = pinus_arizonica[(pinus_arizonica$scientificName=="Pinus arizonica var. cooperi (C. E. Blanco) Farjon" | pinus_arizonica$scientificName=="Pinus cooperi var. ornelasii (Martínez) C. E. Blanco" | pinus_arizonica$scientificName=="Pinus cooperi C. E. Blanco" | pinus_arizonica$scientificName=="Pinus lutea C. E. Blanco ex Martínez"),]

#make the merge
pinus_cooperi_clean = pinus_cooperi[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year")]
pinus_cooperi_arizonica_clean = pinus_cooperi_arizonica[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_cooperi_clean$coordinatePrecision = rep(NA, times=nrow(pinus_cooperi_clean))
pinus_cooperi_clean$coordinateUncertaintyInMeters = rep(NA, times=nrow(pinus_cooperi_clean))
pinus_cooperi_arizonica_clean$coordinatePrecision = rep(NA, times=nrow(pinus_cooperi_arizonica_clean))
pinus_cooperi_final = rbind(pinus_cooperi_clean, pinus_cooperi_arizonica_clean)
str(pinus_cooperi_final)
nrow(pinus_cooperi_final) == nrow(pinus_cooperi_clean) + nrow(pinus_cooperi_arizonica_clean)
unique(pinus_cooperi_final$scientificName)

#write it
write.csv(pinus_cooperi_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_cooperi.csv", row.names = FALSE)

#Pinus fragilissima. It is included in P. taiwanensis. Pinus taiwanensis is similar to P. hwangshanensis, P. luchuensis and P. densiflora and these species are closely related according to phylogenetic analyses using DNA sequence data. Recently, Businský (2003) revisited the morphology of some of these pines and separated some trees in Taiwan as a new species, Pinus fragilissima. Most of the characters evaluated are either similar to those of P. taiwanensis, or they show overlapping states. The seed scales are described as "thin" in the formal description and elsewhere as "fragile" but these are difficult to quantify and may be attributes of the other taxa as well. It also has somewhat longer leaves and only slightly longer seed cones. It is not clear, we will included (buffer). 
pinus_taiwanensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_taiwanensis.csv", header=TRUE)
unique(pinus_taiwanensis$speciesKey) 
unique(pinus_taiwanensis$scientificName) #it is inlcuded fragilissima. here. 
    #Pinus taiwanensis Hayata -> taiwanensis
    #Pinus fragilissima Businský -> fragilissima      
    #Pinus taiwanensis var. taiwanensis -> taiwanensis

pinus_fragilissima = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_fragilissima.csv", header=TRUE)
unique(pinus_fragilissima$speciesKey) 
unique(pinus_fragilissima$scientificName) #Pinus fragilissima Businský
    #Pinus fragilissima Businský

#drop fragilissima from taiwanensis
pinus_taiwanensis_clean = pinus_taiwanensis[!(pinus_taiwanensis$scientificName=="Pinus fragilissima Businský"),]
unique(pinus_taiwanensis_clean$scientificName)

#write it
write.csv(pinus_taiwanensis_clean, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_taiwanensis.csv", row.names = FALSE)

#include these ocurrences in fragilissima
pinus_fragilissima_taiwanensis = pinus_taiwanensis[(pinus_taiwanensis$scientificName=="Pinus fragilissima Businský"),]

#make the merge
pinus_fragilissima_clean = pinus_fragilissima[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year")]
pinus_fragilissima_taiwanensis_clean = pinus_fragilissima_taiwanensis[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_fragilissima_clean$coordinatePrecision = rep(NA, times=nrow(pinus_fragilissima_clean))
pinus_fragilissima_clean$coordinateUncertaintyInMeters = rep(NA, times=nrow(pinus_fragilissima_clean))
pinus_fragilissima_taiwanensis_clean$coordinatePrecision = rep(NA, times=nrow(pinus_fragilissima_taiwanensis_clean))
pinus_fragilissima_final = rbind(pinus_fragilissima_clean, pinus_fragilissima_taiwanensis_clean)
str(pinus_fragilissima_final)
nrow(pinus_fragilissima_final) == nrow(pinus_fragilissima_clean) + nrow(pinus_fragilissima_taiwanensis_clean)
unique(pinus_fragilissima_final$scientificName)

#write it
write.csv(pinus_fragilissima_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_fragilissima.csv", row.names = FALSE)

#Pinus juarezensis. It is included in P.quadrifolia and monophylla. There is a lof of controversy (http://www.conifers.org/pi/Pinus_quadrifolia.php), algunos separan a quadifolia de juarezeneiss, otros dicen que quadrifolia es un hibrido raro no visto casi nunca.. I am not sure, thus I will use the taxa with buffers. P. juarezensis and P. quadrifolia as separated.  
pinus_quadrifolia = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_quadrifolia.csv", header=TRUE)
unique(pinus_quadrifolia$speciesKey) 
unique(pinus_quadrifolia$scientificName) #it is inlcuded juarezensis here. 
    #[1] Pinus quadrifolia Parl. ex Sudw. -> quadrifolia      
    #[2] Pinus juarezensis Lanner -> juarezensis (se poen como sinonimo de quadrifolia porque muchos consideran que lo que se definión como cuadrifolia era un hibrido de monofilla y juarezensis que tenía solo 4 hojas en vez de las 5 de juaarezensis)                 
    #[3] Pinus quadrifolia Parry -> quadrifolia                     
    #[4] Pinus quadrifolia Parry ex Parl. -> quadrigolia             
    #[5] Pinus cembroides var. parryana (Engelm.) Voss -> quadrigolia
    #[6] Pinus parryana Engelm. -> quadrigolia                   

pinus_monophylla = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_monophylla.csv", header=TRUE)
unique(pinus_monophylla$speciesKey) 
unique(pinus_monophylla$scientificName) #It is not included in this species
    #[1] Pinus monophylla Torr. & Frem. -> monophylla                              
    #[2] Pinus edulis var. fallax Little -> monophylla                            
    #[3] Pinus californiarum D.K. Bailey -> monophylla                            
    #[4] Pinus californiarum subsp. fallax (Little) D.K. Bailey -> monophylla      
    #[5] Pinus monophylla var. californiarum (D.K. Bailey) Silba -> monophylla     
    #[6] Pinus cembroides subsp. monophylla (Torr. & Frém.) E. Murray -> monophylla
    #[7] Pinus cembroides subsp. monophylla (Torr. & Frém.) A.E.Murray -> monophylla
    #[8] Pinus cembroides var. monophylla (Torr. & Frém.) Voss -> monophylla     
    #[9] Pinus monophylla Hort. -> monophylla                                       
    #[10] Pinus fallax (Little) Businský -> monophylla

pinus_juarezensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_juarezensis.csv", header=TRUE)
unique(pinus_juarezensis$speciesKey) 
unique(pinus_juarezensis$scientificName) #Pinus juarezensis Lanner
    #Pinus juarezensis Lanner -> juarezensis

#drop juarezensis ocurrences from quadrifola
pinus_quadrifolia_clean = pinus_quadrifolia[!pinus_quadrifolia$scientificName=="Pinus juarezensis Lanner",]
unique(pinus_quadrifolia_clean$scientificName)

write.csv(pinus_quadrifolia_clean, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_quadrifolia.csv", row.names = FALSE)

#include these ocurrences in juarezensis
pinus_juarezensis_quadrifolia = pinus_quadrifolia[pinus_quadrifolia$scientificName=="Pinus juarezensis Lanner",]

#make the merge
pinus_juarezensis_clean = pinus_juarezensis[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year")]
pinus_juarezensis_quadrifolia_clean = pinus_juarezensis_quadrifolia[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_juarezensis_clean$coordinateUncertaintyInMeters = rep(NA, times=nrow(pinus_juarezensis_clean))
pinus_juarezensis_clean$coordinatePrecision = rep(NA, times=nrow(pinus_juarezensis_clean))
pinus_juarezensis_quadrifolia_clean$coordinatePrecision = rep(NA, times=nrow(pinus_juarezensis_quadrifolia_clean))
pinus_juarezensis_final = rbind(pinus_juarezensis_clean, pinus_juarezensis_quadrifolia_clean)
str(pinus_juarezensis_final)
nrow(pinus_juarezensis_final) == nrow(pinus_juarezensis_clean) + nrow(pinus_juarezensis_quadrifolia_clean)
unique(pinus_juarezensis_final$scientificName)

#write it
write.csv(pinus_juarezensis_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_juarezensis.csv", row.names = FALSE)

#Pinus kwangtungensis. It is included in P. fenzeliana and P.morrisonicola. Two varieties. For the type variety, syn: Pinus wangii Hu et W. C. Cheng var. kwangtungensis (Chun et Tsiang) Silba. The second variety is Pinus kwangtungensis var. varifolia Nan Li et Y. C. Zhong, Novon 7: 262. 1997. Many authors do not treat it at the species rank, preferring to lump it with either of the closely related taxa Pinus morrisonicola of Taiwan or Pinus fenzeliana of south China and Vietnam. In this treatment of the group, I have chosen to treat each taxon as a good species, following the line taken by the Flora of China (Wu and Raven 1999) as a recent authoritative treatment. However, I feel the Flora does tend toward "splitting" taxa on the basis of what many authors consider to be minor characters. (http://www.conifers.org/pi/Pinus_kwangtungensis.php). 
pinus_fenzeliana = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_fenzeliana.csv", header=TRUE)
unique(pinus_fenzeliana$speciesKey) 
unique(pinus_fenzeliana$scientificName) #it is inlcuded kwangtungensis here. 
    #[1] Pinus eremitana Businský -> fenzeliana (according gbif and the Catalogue of life)                                
    #[2] Pinus kwangtungensis Chun & Tsiang -> kwangtungensis                      
    #[3] Pinus fenzeliana Hand.-Mazz. -> fenzeliana                             
    #[4] Pinus wangii subsp. kwangtungensis (Chun & Tsiang) Businský -> kwangtungensis
    #[5] Pinus kwangtungensis var. varifolia Nan Li & Y.C.Zhong -> kwangtungensis  
        #The cases indicated as kwangtungensis are considered sinonimous of kwangtungensis, altoguh this name is not accepted for some authors

pinus_morrisonicola = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_morrisonicola.csv", header=TRUE)
unique(pinus_morrisonicola$speciesKey) 
unique(pinus_morrisonicola$scientificName) #it is not included here. 
    #[1] Pinus parviflora var. morrisonicola (Hayata) C.L. Wu -> morrisonicola
    #[2] Pinus morrisonicola Hayata -> morrisonicola                       
    #[3] Pinus formosana Hayata -> morrisonicola                            
    #[4] Pinus uyematsui Hayata -> morrisonicola 

pinus_kwangtungensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_kwangtungensis.csv", header=TRUE)
unique(pinus_kwangtungensis$speciesKey) 
unique(pinus_kwangtungensis$scientificName) #Pinus kwangtungensis Chun & Tsiang
    #Pinus kwangtungensis Chun & Tsiang -> kwangtungensis
    
#drop kwangtungensis from fenzeliana
pinus_fenzeliana_clean = pinus_fenzeliana[!(pinus_fenzeliana$scientificName=="Pinus kwangtungensis Chun & Tsiang" | pinus_fenzeliana$scientificName=="Pinus wangii subsp. kwangtungensis (Chun & Tsiang) Businský" | pinus_fenzeliana$scientificName=="Pinus kwangtungensis var. varifolia Nan Li & Y.C.Zhong"),]
unique(pinus_fenzeliana_clean$scientificName)

#write it 
write.csv(pinus_fenzeliana_clean, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_fenzeliana.csv", row.names = FALSE)

#include these ocurrences in kwangtungensis.
pinus_kwangtungensis_fenzeliana = pinus_fenzeliana[(pinus_fenzeliana$scientificName=="Pinus kwangtungensis Chun & Tsiang" | pinus_fenzeliana$scientificName=="Pinus wangii subsp. kwangtungensis (Chun & Tsiang) Businský" | pinus_fenzeliana$scientificName=="Pinus kwangtungensis var. varifolia Nan Li & Y.C.Zhong"),]

#make the merge
pinus_kwangtungensis_clean = pinus_kwangtungensis[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year")]
pinus_kwangtungensis_fenzeliana_clean = pinus_kwangtungensis_fenzeliana[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_kwangtungensis_clean$coordinatePrecision = rep(NA, times=nrow(pinus_kwangtungensis_clean))
pinus_kwangtungensis_clean$coordinateUncertaintyInMeters = rep(NA, times=nrow(pinus_kwangtungensis_clean))
pinus_kwangtungensis_fenzeliana_clean$coordinatePrecision = rep(NA, times=nrow(pinus_kwangtungensis_fenzeliana_clean))
pinus_kwangtungensis_final = rbind(pinus_kwangtungensis_clean, pinus_kwangtungensis_fenzeliana_clean)
str(pinus_kwangtungensis_final)
nrow(pinus_kwangtungensis_final) == nrow(pinus_kwangtungensis_clean) + nrow(pinus_kwangtungensis_fenzeliana_clean)
unique(pinus_kwangtungensis_final$scientificName)

#write it
write.csv(pinus_kwangtungensis_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_kwangtungensis.csv", row.names = FALSE)

#Pinus maestrensis. It is included in P.occidentalis and P.cubensis. It is not clear. Buffer. 
pinus_cubensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_cubensis.csv", header=TRUE)
unique(pinus_cubensis$speciesKey) 
unique(pinus_cubensis$scientificName) #it is included here. 
    #Pinus cubensis Griseb. -> cubensis
    #Pinus maestrensis Bisse -> maestrensis
    #Pinus wrightii Engelm. -> cubensis
        #The cases indicated as maestrensis are considered sinonimous of maestrensis, altoguh this name is not accepted for some authors

pinus_occidentalis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_occidentalis.csv", header=TRUE)
unique(pinus_occidentalis$speciesKey) 
unique(pinus_occidentalis$scientificName) #it is not included here. 
    #[1] Pinus occidentalis Sw. -> occidentails                 
    #[2] Pinus occidentalis var. baorucoensis Silba -> occidentails


pinus_maestrensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_maestrensis.csv", header=TRUE)
unique(pinus_maestrensis$speciesKey) 
unique(pinus_maestrensis$scientificName) 
    #Pinus maestrensis Bisse -> maestrensis

#drop ocurrences of maestrensis. in pinus_cubensis
pinus_cubensis_clean = pinus_cubensis[!(pinus_cubensis$scientificName=="Pinus maestrensis Bisse"),]
unique(pinus_cubensis_clean$scientificName)

write.csv(pinus_cubensis_clean, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_cubensis.csv", row.names = FALSE)

#include these ocurrences in Pinus maestrensis.
pinus_maestrensis_cubensis = pinus_cubensis[(pinus_cubensis$scientificName=="Pinus maestrensis Bisse"),]

#make the merge
pinus_maestrensis_clean = pinus_maestrensis[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_maestrensis_cubensis_clean = pinus_maestrensis_cubensis[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_maestrensis_clean$coordinatePrecision = rep(NA, times=nrow(pinus_maestrensis_clean))
pinus_maestrensis_cubensis_clean$coordinatePrecision = rep(NA, times=nrow(pinus_maestrensis_cubensis_clean))
pinus_maestrensis_final = rbind(pinus_maestrensis_clean, pinus_maestrensis_cubensis_clean)
str(pinus_maestrensis_final)
nrow(pinus_maestrensis_final) == nrow(pinus_maestrensis_clean) + nrow(pinus_maestrensis_cubensis_clean)
unique(pinus_maestrensis_final$scientificName)

#write it
write.csv(pinus_maestrensis_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_maestrensis.csv", row.names = FALSE)


#Pinus washoensis. It is included in P.ponderosa. There is controversy (http://www.conifers.org/pi/Pinus_washoensis.php). Buffer. 
pinus_ponderosa = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_ponderosa.csv", header=TRUE)
unique(pinus_ponderosa$speciesKey) 
unique(pinus_ponderosa$scientificName) #it is inlcuded washoensis here. 
    #[1] Pinus ponderosa Douglas ex C. Lawson -> ponderosa                                       
    #[2] Pinus ponderosa var. scopulorum Engelm. -> ponderosa                                   
    #[3] Pinus ponderosa subsp. scopulorum (Engelm.) E. Murray -> ponderosa                    
    #[4] Pinus ponderosa Douglas -> ponderosa                                                 
    #[5] Pinus ponderosa P. Lawson & C. Lawson -> ponderosa                                     
    #[6] Pinus ponderosa var. pacifica J. R. Haller & Vivrette -> ponderosa                      
    #[7] Pinus ponderosa var. ponderosa -> ponderosa                                          
    #[8] Pinus ponderosa subsp. scopulorum (Engelm.) A. E. Murray -> ponderosa                  
    #[9] Pinus washoensis Mason & Stockw. -> washoensis                                        
    #[10] Pinus scopulorum (Engelm.) Lemmon -> ponderosa                                          
    #[11] Pinus resinosa Torr. ->  ponderosa                                              
    #[12] Pinus ponderosa var. benthamiana (Hartw.) Vasey ->                             
    #[13] Pinus ponderosa var. washoensis (H. Mason & Stockw.) J. R. Haller & Vivrette -> washoensis
    #[14] Pinus brachyptera Engelm. -> ponderosa                                               
    #[15] Pinus benthamiana Hartw. ->  ponderosa                                                
        #The cases indicated as washoensis are usually considered as synonimun of a intraspecific taxon within ponderosa.

pinus_washoensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_washoensis.csv", header=TRUE)
unique(pinus_washoensis$speciesKey) 
unique(pinus_washoensis$scientificName) #Pinus washoensis Mason & Stockw.

#drop pinus_washoensis from pinus_ponderosa
pinus_ponderosa_clean = pinus_ponderosa[!(pinus_ponderosa$scientificName=="Pinus washoensis Mason & Stockw." | pinus_ponderosa$scientificName=="Pinus ponderosa var. washoensis (H. Mason & Stockw.) J. R. Haller & Vivrette"),]
unique(pinus_ponderosa_clean$scientificName)

write.csv(pinus_ponderosa_clean, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_ponderosa.csv", row.names = FALSE)

#include these ocurrences in pinus_washoensis
pinus_washoensis_ponderosa = pinus_ponderosa[(pinus_ponderosa$scientificName=="Pinus washoensis Mason & Stockw." | pinus_ponderosa$scientificName=="Pinus ponderosa var. washoensis (H. Mason & Stockw.) J. R. Haller & Vivrette"),]

#make the merge
pinus_washoensis_clean = pinus_washoensis[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_washoensis_ponderosa_clean = pinus_washoensis_ponderosa[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinatePrecision", "coordinateUncertaintyInMeters")]
pinus_washoensis_clean$coordinatePrecision = rep(NA, times=nrow(pinus_washoensis_clean))
pinus_washoensis_final = rbind(pinus_washoensis_clean, pinus_washoensis_ponderosa_clean)
str(pinus_washoensis_final)
nrow(pinus_washoensis_final) == nrow(pinus_washoensis_clean) + nrow(pinus_washoensis_ponderosa_clean)
unique(pinus_washoensis_final$scientificName)

#write it
write.csv(pinus_washoensis_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_washoensis.csv", row.names = FALSE)

#Pinus yecorensis. It . There is some doubt about the validity of this taxon. Farjon and Styles (1997) placed it into synonymy with Pinus pseudostrobus var. pseudostrobus and its morphological characters appear to fall fully within the range of variability in that taxon. However, its occurence is more northerly than that of P. pseudostrobus and its ecological setting, in woodlands below the elevation of continuous forest, is somewhat distinct (although our historical maps indicate some overlap). Moreover, a molecular analysis by Gernandt et al. (2009) placed it in a clade with P. douglasiana and P. maximinoi, consistent with the analysis of Felger et al. (2001) who allied it with P. douglasiana on purely morphological grounds. It may prove to be more accurately described as a variety of P. douglasiana but has not yet been named as such.
pinus_pseudostrobus = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_pseudostrobus.csv", header=TRUE)
unique(pinus_pseudostrobus$speciesKey) 
unique(pinus_pseudostrobus$scientificName) #it is inlcuded yecorensis here. 
    #[1] Pinus pseudostrobus Lindl. -> pseudostrobus                        
    #[2] Pinus pseudostrobus var. apulcensis (Lindl.) Shaw -> Pinus apulcensis   
    #[3] Pinus pseudostrobus var. pseudostrobus  -> pseudostrobus           
    #[4] Pinus yecorensis Debreczy & Rácz -> yecorensis                      
    #[5] Pinus oaxacana Mirov (sinonmun of Pinus pseudostrobus var. apulcensis (Lindl.) Shaw) ->  pseudostrobus                           
    #[6] Pinus pseudostrobus var. estevezii Martinez -> pseudostrobus           
    #[7] Pinus pseudostrobus f. protuberans Martinez -> pseudostrobus     
    #[8] Pinus pseudostrobus var. oaxacana (Mirov) S.G. Harrison -> pseudostrobus
    #[9] Pinus estevezii (Martínez) J. P. Perry -> pseudostrobus               
    #[10] Pinus pseudostrobus subsp. oaxacana (Mirov) Silba -> pseudostrobus     
    #[11] Pinus pseudostrobus var. coatepecensis Martinez -> pseudostrobus      
    #[12] Pinus pseudostrobus f. pseudostrobus ->  pseudostrobus                
    #[13] Pinus pseudostrobus subsp. apulcensis (Lindl.) Stead -> pseudostrobus  
    #[14] Pinus nubicola J. P. Perry -> pseudostrobus                          
    #[15] Pinus pseudostrobus var. apulcensis (Lindl.) Martínez -> pseudostrobus 
    #[16] Pinus pseudostrobus var. laubenfelsii Silba  ->  pseudostrobus        
    #[17] Pinus apulcensis Lindl. -> Pinus apulcensis         
        #except yecorensis, the rest of entites seems to avoid the nortwest of Mexico.

pinus_yecorensis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_yecorensis.csv", header=TRUE)
unique(pinus_yecorensis$speciesKey) 
unique(pinus_yecorensis$scientificName) 
    #Pinus yecorensis Debreczy & Rácz -> yecorensis

#drop pinus_yecorensis from pinus_pseudostrobus
pinus_pseudostrobus_clean = pinus_pseudostrobus[!pinus_pseudostrobus$scientificName=="Pinus yecorensis Debreczy & Rácz",]
unique(pinus_pseudostrobus_clean$scientificName)

#write it 
write.csv(pinus_pseudostrobus_clean, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_pseudostrobus.csv", row.names = FALSE)

#include these ocurrences in pinus_yecorensis
pinus_yecorensis_pseudostrobus = pinus_pseudostrobus[pinus_pseudostrobus$scientificName=="Pinus yecorensis Debreczy & Rácz",]

#make the merge
pinus_yecorensis_clean = pinus_yecorensis[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year")]
pinus_yecorensis_pseudostrobus_clean = pinus_yecorensis_pseudostrobus[,c("basisOfRecord", "scientificName", "speciesKey", "lon", "lat", "year", "coordinateUncertaintyInMeters")]
pinus_yecorensis_clean$coordinatePrecision = rep(NA, times=nrow(pinus_yecorensis_clean))
pinus_yecorensis_clean$coordinateUncertaintyInMeters = rep(NA, times=nrow(pinus_yecorensis_clean))
pinus_yecorensis_pseudostrobus_clean$coordinatePrecision = rep(NA, times=nrow(pinus_yecorensis_pseudostrobus_clean))
pinus_yecorensis_final = rbind(pinus_yecorensis_clean, pinus_yecorensis_pseudostrobus_clean)
str(pinus_yecorensis_final)
nrow(pinus_yecorensis_final) == nrow(pinus_yecorensis_clean) + nrow(pinus_yecorensis_pseudostrobus_clean)
unique(pinus_yecorensis_final$scientificName)

#write it
write.csv(pinus_yecorensis_final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_yecorensis.csv", row.names = FALSE)

###Test final species

#extract the name of the files in the folder with final occurrences
list_species_occurrences = list.files("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species", pattern=".csv")

#extract those with error
cases_with_error = NULL
for(i in 1:length(list_species_occurrences)){

    #selected species
    selected_species = list_species_occurrences[i]

    #is error included in the name?
    if(grepl("error", selected_species)){
        cases_with_error = append(cases_with_error, selected_species)
    }
    
}

#check that the total number of species with occurrence data (discounting previous verisons of datasets with errors)
list_species_occurrences = list_species_occurrences[-which(list_species_occurrences %in% cases_with_error)]
length(list_species_occurrences) == 113

#chck that all final species are included in the biancas tree
list_species_occurrences %in% paste(list_species$genus, paste(list_species$specific_epithet, "csv", sep="."), sep="_")

#check that all species of the Bianca's tree are included in the final folder with occurrences
paste(list_species$genus, paste(list_species$specific_epithet, "csv", sep="."), sep="_") %in% list_species_occurrences
#All final species are included in bianca's tree and viceversa

#save the workspace
save.image("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/rdata/taxonomic_discussions.RData")