library(Seurat)
library(gridExtra)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(cowplot)
library(Nebulosa)
library(ggpubr)
library(R.utils)
set.seed(seed = 152727)


#cell labeller models pretrained from scannotater
#default_models <- load_models("default")


#add command line arguments
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])




plotter <- function(s_obj, sample_name,sample_info){
#      ll <- plot_density(s_obj,features = ##c('rna_MSH2','rna_MLH1',"rna_MLH3","rna_MSH3","rna_PMS1",'rna_MSH6','rna_PMS2'),joint = TRUE,combine = #TRUE)
      
sensor2_global <- fread(paste0('../sensor2_results/',sample_name,'_msi_status'))
      
msings_global <- as.numeric(fread(paste0('../msings_results/', sample_name,'_msi_status.MSI_Analysis.txt'))[3,2])*100

pro_global <- fread(paste0('../pro_results/',sample_name,'_msi_status'))
      

plots <-  ggarrange(
  print(FeaturePlot(s_obj,features=c('sensor2'),cols=c('blue','red')) + 
  ggtitle(paste0('MSIsensor2 (',sensor2_global$`%`,' overall)'))),
  
  print(FeaturePlot(s_obj,features=c('msings'),cols=c('blue','red')) + 
  ggtitle(paste0('mSINGS (',msings_global,' overall)'))),
  
   print(FeaturePlot(s_obj,features=c('pro'),cols=c('blue','red')) + 
  ggtitle(paste0('MSIsensor-pro (',pro_global$`%`,' overall)'))),
  
  
  
  
  print(DimPlot(s_obj,group.by = 'numbat',cols=c('blue','red'))+ggtitle(paste0('Numbat'))),
  #print(DimPlot(s_obj)+ggtitle(paste0('Clusters'))),
  print(DimPlot(s_obj,group.by='copykat',cols=c('red','blue','grey'))+ggtitle(paste0('Copykat')))
  )

print(annotate_figure(plots,top=text_grob(paste0(sample_info$filename,'; ',sample_info$site,'; ',
                                                sample_info$msi_status,'; ',sample_info$msi_test),
color='black',face='bold',size=20)))
#print(ggarrange(FeaturePlot(s_obj,features=c('rna_MSH2')),
#                FeaturePlot(s_obj,features='rna_MSH3'),
#                  FeaturePlot(s_obj,features=c('rna_MSH6')),
#                  FeaturePlot(s_obj,features=c('rna_MLH1')),
#                FeaturePlot(s_obj,features='rna_MLH3'),
##                FeaturePlot(s_obj,features=c('rna_PMS2')),
#                FeaturePlot(s_obj,features='rna_KRAS'),
#                FeaturePlot(s_obj,features='rna_TP53'),
#                FeaturePlot(s_obj,features='rna_BRAF')
#                ))
                #  common.legend = TRUE,
                 # legend = 'right'))
#print(ggarrange(FeaturePlot(s_obj,features='MSH2'),
##FeaturePlot(s_obj,features='MSH3'),
#FeaturePlot(s_obj,features='MSH6'),
##FeaturePlot(s_obj,features='MLH1'),
#FeaturePlot(s_obj,features='MLH3')))
        ##common.legend=TRUE,
        #legend='right'))
#print(ggarrange(FeaturePlot(s_obj,features='PMS1'),
#FeaturePlot(s_obj,features='PMS2'),
#FeaturePlot(s_obj,features='KRAS'),
#FeaturePlot(s_obj,features='TP53'),
#FeaturePlot(s_obj,features='BRAF')))
        #common.legend=TRUE,
        #legend='right'))


#print(ll)
  
}

#to print violin plot VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
plotter2 <- function(s_obj, sample_name){
genes <- colnames(s_obj[[]])
genes <- genes[12:length(genes)]

plots <- FeaturePlot(s_obj,features=genes,combine=FALSE,cols=c('khaki','red'))



grobs <-marrangeGrob(grobs = plots, nrow=3, ncol=3)
return(grobs)

#){
#  print(FeaturePlot(s_obj, features = gene))
}
  
#}

key <- fread('../manifests/final_key.tsv')
samp_list <- fread('../manifests/all_samples.tsv',header=FALSE)
#or better yet
samp_list <- list.files(path='../annotated_h5/',pattern = '.rds',full.names=FALSE)
samp_list <- gsub(pattern = '.rds',replacement = '',x = samp_list)

#sanity check with just one
#samp_list <- samp_list[1]

sample <- sample_name
  
sample_info <- filter(key, filename == sample )
if(nrow(sample_info) == 0){
  sample_info <- data.frame(patient_id = sample,
                            site ='tumor',
                            filename = sample,
                            msi_status = 'MSI-H',
                            msi_test='IHC/PCR')

  
}
assign(x = paste0(sample),value = readRDS(paste0('../annotated_h5/',sample,'.rds')))

ggsave(paste0('../images/',sample,'_msi_plots.pdf'),plotter(eval(readRDS(paste0('../annotated_h5/',sample,'.rds'))),
             paste0(sample),
              sample_info=sample_info),width=3, height=3, units="in", scale=3)

ggsave(paste0('../images/',sample,'_mut_plots.pdf'),plotter2(eval(readRDS(paste0('../annotated_h5/',sample,'.rds'))),
          paste0(sample)),width=3, height=3, units="in", scale=3)

