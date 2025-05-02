#first iteration based on function found here: https://bioinformatics-core-shared-training.github.io/SingleCell_RNASeq_May23/UnivCambridge_ScRnaSeqIntro_Base/Markdowns/101-seurat_part2.html
#also see the merge function from the seurat v4.3 archive https://satijalab.org/seurat/archive/v4.3/merge

#set seed (remove this for random sampling)
#set.seed(9)

#load libs
library(Seurat)
library(subSeq)
library(data.table)
library(tidyverse)
library(R.utils)
library(DropletUtils)


#add command line arguments
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
manifest <- as.character(argus[1])
sample1_name <- as.character(argus[2])
sample2_name <- as.character(argus[3])
runs <- as.numeric(argus[4])


#load data 
print(paste0("using manifest: ","../manifests/",manifest))

key <- fread(paste0('../manifests/',manifest))


#define funcs
cell_sampler <- function(sample_name,msi_status,prop_cancer,prop_normal){

#description of paramters
#sample_name; give sample name of s_obj to be used
#prop_cancer: give proportion of cancer cells to be kept ex: .9
#prop_normal give proportion of normal cells to be kept

#read in annotated seurat object. Must have scATOMIC annotations
s_obj <- readRDS(paste0('../annotated_h5/',sample_name,'.rds'))

#add gt for downstream cell tracking
s_obj$ground_truth <- paste(msi_status,s_obj$pan_cancer_cluster,sep="_")

cancer_cells <- subset(s_obj,subset=pan_cancer_cluster == "Cancer")
normal_cells <- subset(s_obj,subset=pan_cancer_cluster == "Normal")

num_cancer_cells <- as.numeric(ncol(cancer_cells) * prop_cancer)

num_normal_cells <- as.numeric(ncol(normal_cells) * prop_normal)

tot_cells = as.numeric(num_cancer_cells + num_normal_cells)

#currently not operation but more accurate classification of "normal" or "cancer" where they are defined to be from a "cancer cluster"
#These annotations are only currently generated as part of the integrated sample portion of the pipeline.
##num_norm_cells <- ncol(subset(s_obj,subset=pan_cancer_cluster == "Normal" & percent_cancer < 70))
##num_cancer_cells <- ncol(subset(s_obj,subset=pan_cancer_cluster == "Cancer" & percent_cancer >= 70))



#subsample seurat object based on desired number of cancer cells
sampled_cancer <- cancer_cells[, sample(colnames(cancer_cells), size = num_cancer_cells, replace=F)]
sampled_normal <- normal_cells[, sample(colnames(normal_cells), size = num_normal_cells, replace=F)]

#merge subsamples
#can add cell idents with add.cell.ids
subbed <- merge(sampled_cancer,y=sampled_normal)
subbed[["RNA"]] <- JoinLayers(subbed[["RNA"]])

#add sampling information to s_obj metadta
subbed$tot_cells <- tot_cells
subbed$num_cancer_cells <- num_cancer_cells
subbed$num_normal_cells <- num_normal_cells
subbed$msi_grounded <- msi_status

return(subbed)
}


#define mix intervals 
mix_intervals <- seq(from = 0.1, to = .9, by = .1)

#perform multiple runs (to be averaged at the end. Adjust the cv variable to change the number of runs)
all_mixing_tables <- data.frame()

for(cv in 1:runs){
print(paste0('starting run ', cv))

#setting seed for each run to the n'th run. Easy way to keep track of seeds and things. Probably should switch to random seed generator though...
set.seed(cv)






#create mixing table to generate all possible combinations of MSIH and MSS cells
mixing_table <- data.frame(prop_msih_cancer=mix_intervals,prop_mss_cancer=rev(mix_intervals),
			   prop_msih_norm=mix_intervals,prop_mss_norm=rev(mix_intervals))





#generate msih and mss lists of seurat objects with varying mixes of cells based on mixing table.

print(paste("creating objects for sample",sample1_name))
msih_objs <- apply(X=mixing_table,function(x) cell_sampler(sample_name=sample1_name,msi_status='msih',x[1],x[3]),MARGIN=1)

print(paste("creating objects for sample",sample2_name))
mss_objs <- apply(X=mixing_table,function(x) cell_sampler(sample_name=sample2_name,msi_status='mss',x[2],x[4]),MARGIN=1)


#add ground truth metadata
for(i in 1:nrow(mixing_table)){
mixing_table$mss_num_norm[i] <-floor(as.numeric(mss_objs[[i]]$num_normal_cells[1]))
mixing_table$mss_num_cancer[i] <-floor(as.numeric(mss_objs[[i]]$num_cancer_cells[1]))
mixing_table$mss_tot_cells[i] <- floor(as.numeric(mss_objs[[i]]$tot_cells[1]))
mixing_table$msih_num_norm[i] <- floor(as.numeric(msih_objs[[i]]$num_normal_cells[1]))
mixing_table$msih_num_cancer[i] <- floor(as.numeric(msih_objs[[i]]$num_cancer_cells[1]))
mixing_table$msih_tot_cells[i] <-floor(as.numeric(msih_objs[[i]]$tot_cells[1]))
}



#merge mixes (so .9 msih mixes with .1 mss)
print("mixing samples together")
all_tomorrows <- Map(merge,msih_objs,mss_objs)
for(i in 1:length(all_tomorrows)){
all_tomorrows[[i]][["RNA"]] <- JoinLayers(all_tomorrows[[i]][["RNA"]])
}

#write out artifical_manifest (if we want to do this earlier in the script, we'll have to change the apply function
#as it coerces the dataframe into a single variable type
mixing_table$sample_id <-c(paste0("mix_",1:9,'_',cv))
mixing_table$file_prefix <- mixing_table$sample_id

fwrite(x=mixing_table,file=paste0('../manifests/mixing_table_manifest_',cv,'.tsv'),sep='\t')

all_mixing_tables <- rbind(all_mixing_tables,mixing_table)


#create artificial 10x samples for use with SC-MSI pipeline
for(i in 1:length(all_tomorrows)){
write10xCounts(path=paste0('../artificial_samples/mix',i,'_',cv),
x=all_tomorrows[[i]]@assays$RNA@layers$counts,
version="3",barcodes = colnames(all_tomorrows[[i]]),
  gene.id = rownames(all_tomorrows[[i]]),
  gene.symbol = rownames(all_tomorrows[[i]]),
  gene.type="Gene Expression", overwrite=TRUE, type="auto")

file.copy(from=paste0('../artificial_samples/mix',i,'_',cv,'/barcodes.tsv.gz'),
	  to=paste0('../gsm_samps/mix_',i,'_',cv,'_barcodes.tsv.gz'),overwrite=TRUE)
file.copy(from=paste0('../artificial_samples/mix',i,'_',cv,'/features.tsv.gz'),
          to=paste0('../gsm_samps/mix_',i,'_',cv,'_features.tsv.gz'),overwrite=TRUE)
file.copy(from=paste0('../artificial_samples/mix',i,'_',cv,'/matrix.mtx.gz'),
          to=paste0('../gsm_samps/mix_',i,'_',cv,'_matrix.mtx.gz'),overwrite=TRUE)

#save ground truth annotations
print('saving ground truth objects')
saveRDS(object=all_tomorrows[[i]],file=paste0('../temp/mix_',i,'_',cv,'ground_truth.rds'))




}

print(paste0("Finished Run ",cv))
}

#write out artifical_manifest containing all artificial samples from all runs 
fwrite(x=all_mixing_tables,file=paste0('../manifests/all_mixing_tables_manifest.tsv'),sep='\t')


#write out ALL samples names
writeLines(text=all_mixing_tables$sample_id,paste0("../manifests/all_artificial_samples.tsv"),sep='\n')
#copy sample names to create patient list (redundant, so might want to remove this when I have time
file.copy(from=paste0("../manifests/all_artificial_samples.tsv"),
          to=paste0("../manifests/all_mix_patients.txt"),overwrite=TRUE)




#also create a new mixing_key for Snakemake
mix_key_template <- fread('../manifests/mix_key.tsv')
mix_key_rep <- do.call("rbind", replicate(runs, mix_key_template, simplify = FALSE))

mix_key_rep$sample_id <- paste0(mix_key_rep$sample_id,'_',1:runs)
mix_key_rep$file_prefix <- paste0(mix_key_rep$file_prefix,'_',1:runs)
mix_key_rep$patient_id <- mix_key_rep$sample_id

fwrite(mix_key_rep,paste0('../manifests/mix_key_cv.tsv'),sep='\t')







