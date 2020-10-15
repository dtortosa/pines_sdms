#set working directory
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/with_proportions")

#list species
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)
str(list_species)
summary(list_species)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false

#drop discolor (problem tazonomy, no diferetiaced from cembriodes). Also remove tecunumi y jaslicana, for which we didn't create maps in the beginning
epithet_species_list=epithet_species_list[which(!epithet_species_list %in% c("discolor", "jaliscana", "tecunumanii"))]

#check it
!c("discolor", "jaliscana", "tecunumanii") %in% epithet_species_list
length(epithet_species_list) == 112

#create a vector with the complete name of the pdf
final_list_species = paste(epithet_species_list, "_with_proportions.pdf", sep="")

#copy a necesary program 
#system("cp /Users/diegosalazar/cpdf-binaries-master/OSX-Intel/cpdf /usr/bin/cpdf")

#add this package to .bash_profile to use it in any directory
#~/.bash_profile
#export PATH=$PATH:/usr/bin/cpdf

#create sequence of numbers to combine each pair of species into one figure (speices1 and species2 into figure 1, species3 and speceis4 into figure 2). Therefore, for each two species we select one number, until 112 (1,3,5,...)
seq_binding_species = seq(1,length(final_list_species),2)
length(seq_binding_species) == 56 #56 pages
summary(seq_binding_species%%2==0) #no even number, only 1,3,5,7,9,11...

#loop for binding each pair of species
for(i in 1:length(seq_binding_species)){

    #select the index for selecting species to bind in the pdf
    index_selecting_species = seq_binding_species[i]

    #selected species
    first_species = final_list_species[index_selecting_species]
    second_species = final_list_species[index_selecting_species+1]

    #set the path for saving images
    output_path = "/Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/with_proportions_2_species_per_page"
    #if the number figure for the pdf file
    number_figure = i
 
    #set the number of page for the PhD
    number_pages=253+i

    #load the
    system(paste("cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/with_proportions;", #set working directory
        "pdfjam ", first_species, " ", second_species, " --paper a4paper --nup 1x2 --no-landscape --scale 0.925 --trim '0cm 8.5cm 0cm 0cm' --offset '0cm 0.7cm' --delta '0cm -0.1cm' --outfile ", output_path, "/figure_", number_figure, ".pdf;", #take the pdf of first and second species, bind them into a4paper (--paper a4paper), in only one column (--nup 1x2 --no-landscape), slightly reduce the size of both figures (--scale 0.925; scale = 1 would be the same size), remove blanck space from the bottom (--trim '0cm 8.5cm 0cm 0cm'; left, bottom, right and top), put more separaton between the bottom-top margins and the figures (--offset '0cm 0.7cm'; offX, offy; see "https://osl.ugr.es/CTAN/macros/latex/contrib/pdfpages/pdfpages.pdf"), control the separation btween figures (--delta '0cm -0.1cm'; see the previous link), and save in the output path using the number figure (from 1 to 56; 1:length(seq_binding_species).
        "cpdf -add-text '", number_pages, "' -pos-center '504 38.5' -font 'Times-Roman' -font-size 11.25 ", output_path, "/figure_", number_figure, ".pdf", " -o ",  output_path, "/figure_", number_figure, ".pdf", sep="")) #add number pages (-add-text) in the bottomright position (-pos-center '501 35'; first number ix the X, the second is the Y), in Timres New Roman (-font 'Times-Roman'), select the figure obtained from pdfjm and overwritte it. See "/Users/diegosalazar/cpdf-binaries-master/cpdfmanual.pdf".
}

#check the number is 56
list_pdfs = list.files("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/with_proportions_2_species_per_page/", pattern=".pdf", full.names = TRUE)

#reorder files according to numerically order
library(gtools)
list_pdfs_order = mixedsort(list_pdfs)
length(list_pdfs_order)==length(seq_binding_species) #TRUE

#bind all pdfs into one single pdf
require(staplr)
staple_pdf(input_files=list_pdfs_order, output_filepath = "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/full_2_species_per_page/full_with_proportions.pdf")

#copy the resulting file (full_with_proportions.pdf) into the directory for binding PhD chapters
system("cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/full_2_species_per_page; cp full_with_proportions.pdf /Volumes/GoogleDrive/My\\ Drive/science/phd/capitulos_maquetados/pdfs_to_bind_v2/supplementary_chapter_3_part_2_v2.pdf")