###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#library
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel


#######################################
###########Load variables##############
####################################### 
list_variables = list.files("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals", full.names=TRUE, pattern=".asc")
require(raster)
variables = stack(list_variables)


##################################
########Extract D2################
##################################

#requireed functions
Dsquared = function(model, adjust = TRUE) {
  # version 1.1 (13 Aug 2013)
  # calculates the explained deviance of a GLM
  # model: a model object of class "glm"
  # adjust: logical, whether or not to use the adjusted deviance taking into acount the nr of observations and parameters (Weisberg 1980; Guisan & Zimmermann 2000). More data used for fitting enhance d2, whilst more parameters lower it.
  d2 <- (model$null.deviance - model$deviance) / model$null.deviance
  if (adjust) {
    n <- length(model$fitted.values)
    p <- length(model$coefficients)
    d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2)
  }
  return(d2)
} #D2. If you have problem with this you can use Dsquared function of modEvA (http://modeva.r-forge.r-project.org/). 

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

#make a function
d2_total = function(i){
    print(i)
    ocurrences_data = read.csv(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(i, "complete.presences.csv", sep="_"), sep="/"), header=TRUE)

    #extract values of variables in each point for creating predictors variables
    predictors = as.data.frame(extract(x=variables, y=ocurrences_data[,c("longitude", "latitude")]))

    #create a unique data frame with response (points) and predictors (environmental variables)
    data = cbind(ocurrences_data, predictors)

    #loop for extract D2 for each variable 
    d2=NULL
    for (k in 1:ncol(predictors)){
        model = glm(data$presence ~ poly(data[,c(5:30)][,k], 2), family=quasibinomial(link=logit), weights=data$precision_weight) #I have used quasibinomial becasue I have weighted the response variable with precision_weight variable, thus we have decimals between 0 and 1. With binomial family we would obtain the following error "fitted probabilities numerically 0 or 1 occurred", but the results are the same.
        d2 = append(d2, Dsquared(model))
    }

    d2_data_frame = data.frame()
    d2_data_frame = rbind(d2_data_frame, d2)
    d2_data_frame = cbind(i, d2_data_frame)
    names(d2_data_frame) = c("species", names(variables))
    d2_data_frame
}

#create a vector with all species
species = epithet_species_list

# set up cluster
clust <- makeCluster(4) 
registerDoParallel(clust)

# run for all species
d2.final = foreach(i = species, .packages = "raster", .combine="rbind") %dopar% { 
    d2_total(i = i)
} #the "stringr" package is used for "str_split_fixed" function

#stop the cluster 
stopCluster(clust)

#check if all is correct
str(d2.final)

#write the rank
write.csv(d2.final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/clustering_results/d2_species_variables.csv", row.names=FALSE)

#create the rank 
d2_rank  = data.frame(t(apply(-d2.final[-1], 1, rank, ties.method='min')))
d2_rank = cbind(d2.final$species, d2_rank)
names(d2_rank)[which(names(d2_rank)=="d2.final$species")] <- "species" 
write.csv(d2_rank, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/clustering_results/d2_rank_species_variables.csv", row.names=FALSE)

### comparison con results run during WSL stay
#d2.final
d2.final_ancient = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/clustering_results/d2_species_variables_wsl.csv", sep=",", header=T) #load d2 values per variable calculated at WSL
identical(d2.final, d2.final_ancient) #check identical. No TRUE because of the decimals 
sum_differences = NULL #calculate the difference between both datasets
for(i in 1:nrow(d2.final)){ #for each row of d2.final

    #select the [i] row without species column for both data.sets
    selected_row_new = d2.final[i,-1]
    selected_row_ancient = d2.final_ancient[i,-1]

    #calculate differences between both data set for each cell of the [i] row. Round to 5 and sum. There are liiiitle differences, 1.5*10^-6 y so on... 
    differences = sum(round(selected_row_new - selected_row_ancient,5))

    #save it
    sum_differences = append(sum_differences, differences)
}
length(sum_differences) #113, one per species
unique(sum_differences) #all zeros

#d.2_rank
d2_rank_ancient = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/clustering_results/d2_rank_species_variables_wsl.csv", sep=",", header=T) #load d2.rank per species calculated at WSL
identical(d2_rank, d2_rank_ancient) #identical. In this case, we don't have decimals


#############################
######Cluster analysys#######
#############################

#load the data set with the rank
d2_rank = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/clustering_results/d2_rank_species_variables.csv", header=TRUE, row.names=NULL)
str(d2_rank)

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

#Load the grid where all the species of the group are stored.
require(raster)
distribution = NULL
for (i in epithet_species_list){
    distri = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(i, "01.img", sep="_"), sep="_"))
    distribution = append(distribution, distri)
}

all_species_grid = stack(distribution)
nlayers(all_species_grid) #113, on for each species
names(all_species_grid) <- epithet_species_list #change the name for match with name sof rank data set speices

#calculate euclidean distance between species 
dendogram_table = d2_rank[-1] #select all columns except species column
rownames(dendogram_table) = epithet_species_list #put species names as row names
species.dist<-dist(abs(dendogram_table), method="euclidean") #We will use the typical euclidean distance: The sum of the squared differences for each value of two vectors ((x[1]-y[1])+(x[2]-y[2])....). Usual distance between the two vectors (2 norm aka L_2), sqrt(sum((x_i - y_i)^2)).

#make a dendrongram with hclust
species_dendogram<-hclust(species.dist) 
plot(species_dendogram)

#save it
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/clustering_results/species_dendogram.png", width=1500, height=800, pointsize=30 )
plot(species_dendogram, cex=0.5)
dev.off()

#load the package neccesary
require(NbClust)
source("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/nbclust_modified.R") #load the nbclust function modified. The rationale is that I obtain an error with complete method: "Error in solve.default(W) : system is computationally singular: reciprocal condition number = 1.12144e-16". 
#solve(w) is used internally by nbclust to calculate the index of friedman, one of the index for evaluate the quality of the number of clusters. Thus is not ver very important because we have other index. The problem is caused because solve is a numerical way to calculate the inverse of the covariance matrix. Unfortunately, if some of the numbers used in the inverse calculation are very small, it assumes that they are zero, leading to the assumption that it is a singular matrix. This is why it specifies that they are computationally singular, because the matrix might not be singular given a different tolerance. The solution is use a smaller tolerance in solve functuin, like solve(..., tol = 1e-19). This should be fine since you get reciprocal condition number = 1.16873e-16 (see http://stackoverflow.com/questions/21451664/system-is-computationally-singular-error). 


#create a function for make clusters and plot the results
species_clustering = function(min_groups, max_groups, method){

    #load the data set with the rank
    d2_rank = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/clustering_results/d2_rank_species_variables.csv", header=TRUE, row.names=NULL)
    str(d2_rank)

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

    #Load the grid where all the species of the group are stored.
    require(raster)
    distribution = NULL
    for (i in epithet_species_list){
        distri = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(i, "01.img", sep="_"), sep="_"))
        distribution = append(distribution, distri)
    }

    all_species_grid = stack(distribution)
    nlayers(all_species_grid) #113, on for each species
    names(all_species_grid) <- epithet_species_list #change the name for match with name sof rank data set speices

    # Create the empty results list
    final_results <- list()

    #run the cluster
    cluster = NbClust(data=d2_rank[-1], #data: matrix or dataset.
        diss=NULL, #dissimilarity matrix to be used (matriz de distancias). You can used one make it before, or if NULL the matrix will be created by means of the distance methos indicated in the next argument. I have test if there is differences between make the matrix before or use the matrix created by nbclust with kmeans method, and there is not differences (see kmeans_2_10_g_3_distances_nbclust.png). 
        distance = "euclidean", #distance: the distance measure to be used to compute the dissimilarity matrix (euclidean, maximum, etc..). It is not neccesary if we introduce the matrix created before. 
        min.nc=min_groups, #minimal number of clusters, between 1 and (number of objects - 1)
        max.nc=max_groups, #maximal number of clusters, between 2 and (number of objects - 1), greater or equal to min.nc. By default, max.nc=15.
        method = method, #the cluster analysis method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans". Achilleas told me thay I should use complete (The distance between two clusters and is the maximum distance between two points x and y). Info of kmeans (see documentation of NbClust function and diapo7 in http://slideplayer.com/slide/5178523/).
        index = "all") #the index to be calculated. They inform about the quality of each proposed clustering (2,3,4,... clusters). If you use more than 1 index, the function summarize how many index consider each number of cluster as optimal, thus how much more index use better). index="all" include a lot of index, but not all. For more info see "https://mail.google.com/mail/u/0/#starred/162867b038d6aa92". 

    #Extracting clustering parameters
    cluster$All.index #results of all index for each number of clusters. Values of indices for each partition of the dataset obtained with a number of clusters between min.nc and max.nc.
    cluster$Best.nc #Best number of clusters proposed by each index and the corresponding index value.
    cluster$All.CriticalValues #Critical values of some indices for each partition obtained with a number of clusters between min.nc and max.nc.

    #Extracting clusterings GROUP details
    groups <- cluster$Best.partition #Partition that corresponds to the best number of clusters
    print(table(groups)) #number of species by group
    num_groups <- length(unique(groups)) #numbe of groups
    group_ids <- unique(groups) #id of the groups
    final_results[[1]] <- cluster #save all results of cluster in final_results
    d2_rank = cbind(d2_rank, groups) #combine d2_rank with the group factor

    #export the group variable with the species names
    write.csv(d2_rank, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_tables/", method, "_", min_groups ,"_", max_groups,"_","g","_", num_groups,".csv",sep=""), row.names = FALSE)

    # Creating the spatial distribution of groups
    groups_final <- list() 
    species_grids <- list()
    for(i in 1:length(group_ids)){ #for each group
        group_grids <- list()
        groups_final[[i]] <- d2_rank[which(groups == group_ids[i]),] #select rows of d2_rank of the group [i]
        sp_names <- as.character(d2_rank[which(groups == group_ids[i]),]$species) #select the speceis of the group [i]
        species_grids[[i]] <- all_species_grid[[match(sp_names,names(all_species_grid))]] #select the distribution grids of the species of the group [i]
    }
    final_results[[2]] <- species_grids #species_grids is a list with the distribution rasters of each species

    # Exporting the rasters with the spatial distribution of groups
    for(j in 1:num_groups){
        print("Exporting raster")
        print(j)
        print("**************")
        writeRaster(species_grids[[j]], paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_grids/", "_", min_groups ,"_", max_groups,"_","g","_", num_groups,"_group_",j ,"_species_",nlayers(species_grids[[j]]), method,".grd",sep=""), overwrite=T)
    } #for each group plot all distribution rasters

    # Plotting the spatial distribution of groups
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/figures/", method, "_", min_groups ,"_", max_groups,"_","g","_", num_groups, ".png",sep=""),width= 3000, height = 2000, bg="white",res=150)
    par(mfcol=c(2,ceiling(num_groups/2)))
    for(j in 1:num_groups){ #for each group
        print("**************")
        print(j)
        print("**************")
        if(length(species_grids[[j]])!=1){ #if there more than 1 layer
          #writeRaster(species_grids[[j]],paste("S:/psomas/requests/scenario_comparison/species_clustering/_grids/",output_name, "_", min_groups ,"_", max_groups,"_","g","_",num_groups,"_group_",j,"_species_",nlayers(species_grids[[j]]),".grd",sep=""))
          plot(sum(species_grids[[j]]>0),main="") #sum all raster of each group, but only the areas with value higher than 1. Plot it. 
          title(main=paste("Species:",nlayers(species_grids[[j]]), sep=" "), cex.main=3)          
        }
        if(length(species_grids[[j]])==1){ #if there is only one layer
          #writeRaster(species_grids[[j]],paste("S:/psomas/requests/scenario_comparison/species_clustering/_grids/",output_name, "_", min_groups ,"_", max_groups,"_","g","_",num_groups,"_group_",j,"_species_",nlayers(species_grids[[j]]),".grd",sep=""))
          plot(species_grids[[j]]>0,main=paste("Species:",nlayers(species_grids[[j]]),sep=" "))
        }
    }
    dev.off()
}

###Paralellize
# set up cluster
clust <- makeCluster(1) #we use fork, because we don't want to copy for each time the rank and the raster of the species. 
registerDoParallel(clust)

#methods to make the clustering
methods = c("average", "centroid", "complete", "kmeans", "mcquitty", "median", "single", "ward.D", "ward.D2")

#run for all species
foreach(method = methods) %dopar% { 
    species_clustering(min_groups=2, max_groups=10, method=method)
} #the "stringr" package is used for "str_split_fixed" function

#stop the cluster 
stopCluster(clust)