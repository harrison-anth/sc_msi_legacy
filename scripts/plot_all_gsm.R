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

plots <-  marrangeGrob(list(print(FeaturePlot(s_obj,features=c('sensor_rna_prob'),cols=c('blue','red')) +
 ggtitle('sensor_rna probs')),
  
print(DimPlot(s_obj,group.by='pan_cancer_cluster',cols=c('red','blue','grey'))+ggtitle('ATOMIC')),
print(DimPlot(s_obj,group.by='sensor_rna_status',cols=c('red','blue'))+ggtitle('sensor_rna'))),
nrow=3,ncol=2)


#plots_final <- print(annotate_figure(plots,top=text_grob(paste0(sample_info$filename,'; ',sample_info$site,'; ',
#                                                sample_info$msi_status,'; ',sample_info$msi_test),
#color='black',face='bold',size=20)))

double_plot <- ggarrange(print(DimPlot(s_obj)+ggtitle(paste0('Clusters'))),
print(DimPlot(s_obj, group.by = "scATOMIC_pred") + ggtitle("ATOMIC cells")+
labs(fill="cell types")+
guides(colour = guide_legend(override.aes = list(size=2.5),ncol=1)))+ theme(legend.text=element_text(size=8)),nrow=1,ncol=1)

double_double <- marrangeGrob(double_plot,ncol=1,nrow=1)

ggsave(paste0('../images/',sample,'_cell_types.pdf'),double_double)

ggsave(paste0('../images/',sample,'_msi_plots.pdf'),plots,width=3, height=3, units="in", scale=3)


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
  
plotter2 <- function(s_obj,s_obj2,s_obj3, sample_name){
genes <- colnames(s_obj3[[]])
genes <- genes[12:length(genes)]

#depth
plots <- FeaturePlot(s_obj,features=genes,combine=FALSE,cols=c('khaki','red'))
#+
#plot_annotation(
#  subtitle = 'Total read depth'
#) 
#af
plots2 <- FeaturePlot(s_obj2,features=genes,combine=FALSE,cols=c('khaki','red'))
#+
#plot_annotation(
#  subtitle = 'Mutant allele frequency'
#)

#alt
plots3 <- FeaturePlot(s_obj3,features=genes,combine=FALSE,cols=c('khaki','red'))
#+
#plot_annotation(
#  subtitle = 'Mutant allele count'
#)


interleaved_plots <- list()
for (i in 1:length(plots)/3){
interleaved_plots <- append(interleaved_plots, plots2[i:(i+2)])
interleaved_plots <- append(interleaved_plots, plots3[i:(i+2)])
interleaved_plots <- append(interleaved_plots, plots[i:(i+2)])
}




grobs <- marrangeGrob(grobs = interleaved_plots, nrow=3, ncol=3)

ggsave(paste0('../images/',sample,'_mut_plots.pdf'),grobs,width=3, height=3, units="in", scale=3)

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





#depth <- assign(x = paste0(sample),value = readRDS(paste0('../annotated_h5/',sample,'.rds')))
plotter(eval(readRDS(paste0('../annotated_h5/',sample,'.rds'))),
             paste0(sample),
              sample_info=sample_info)

#plotter2(eval(readRDS(paste0('../annotated_h5/',sample,'.rds'))),
#eval(readRDS(paste0('../annotated_h5/',sample,'_af.rds'))),
#eval(readRDS(paste0('../annotated_h5/',sample,'_alt_count.rds'))),
#paste0(sample))

