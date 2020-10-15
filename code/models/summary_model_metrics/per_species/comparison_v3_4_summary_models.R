#########################################################################################
######### COMPARING SUMMARY STATISTICS BETWEEN VERSION 3 AND 4 OF THE MODELS ############
#########################################################################################

#load the two versions of the summary statistics
summary_statistics_v3 = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/evaluations/medians_evaluations_v3.csv", header=T, sep=",")
summary_statistics_v4 = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/evaluations/medians_evaluations_v4.csv", header=T, sep=",")

#check that colnames are similar for sd only
colnames(summary_statistics_v3) == colnames(summary_statistics_v4)

#check that the column number is the same
ncol(summary_statistics_v3) == ncol(summary_statistics_v4)

#plot median of third version against mean of the fourth version
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/evaluations/comparsion_summary_models_v3_v4.pdf")
par(mfcol=c(4,2), mai=c(0.3,0.3,0.3,0.3))
for(i in 2:ncol(summary_statistics_v3)){ #for each column with statistics

	#selected [i] columns for the two version of statistics
	summary_statistics_v3_selected_column = summary_statistics_v3[,i]
	summary_statistics_v4_selected_column = summary_statistics_v4[,i]

	#extract the column names
	selected_model_statistic = colnames(summary_statistics_v3)[i]

	#extract the statistic name, which will be removed from the final name
	selected_statistic = strsplit(selected_model_statistic, split="_")[[1]][3]

	#modify the final title
	if(selected_statistic == "sd"){ #if the statistic is sd
		#do not change anything 
		selected_model_statistic_final = selected_model_statistic
	} else { #if not and hence we have mean and median

		#remove the statistic name
		selected_model_statistic_final = strsplit(selected_model_statistic, split=selected_statistic)

		#add to the final name "median_vs_mean"
		selected_model_statistic_final = paste(selected_model_statistic_final, "median_vs_mean", sep="")
	}

	#plot both versions
	plot(x=summary_statistics_v3_selected_column, y=summary_statistics_v4_selected_column, main=selected_model_statistic_final)
	
	#print the selected model and the spearman correlation test
	print("#######################")
	print(selected_model_statistic_final)
	print("#######################")	
	print(cor.test(summary_statistics_v3_selected_column, summary_statistics_v4_selected_column, method="spearman"))
}
dev.off() #SD IS EXACTLY THE SAME, MEAN AND MEDIAN ARE VERY CORRELATED.  


#save workspace
save.image("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/evaluations/comparsion_summary_models_v3_v4.RData")