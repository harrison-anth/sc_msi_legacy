library(pagoda2)
library(numbat)
setwd('C:/Users/Harrison Anthony/Desktop/git_repos/single_msi/markdown_files/')
cd <- read.10x.matrices('../cell_ranger_output/SRR12387937/outs/ell_ranger_output/SRR12387937/outs/raw_feature_bc_matrix/')
counts <- gene.vs.molecule.cell.filter(cd,min.cell.size=500,plot=FALSE)
rownames(counts) <- make.unique(rownames(counts))
counts <- counts[rowSums(counts)>=10,]
counts@Dimnames[[2]] <- gsub(x=counts@Dimnames[[2]],pattern = 'one_',replacement = '')



r <- Pagoda2$new(counts,log.scale=TRUE)
r$adjustVariance(plot=TRUE, gam.k=10)
r$calculatePcaReduction(nPcs=50, n.odgenes=3e3)
r$makeKnnGraph(k=40, type='PCA', center=TRUE, distance='cosine')
r$getKnnClusters(method=infomap.community, type='PCA')
M <- 30
r$getEmbedding(type='PCA', embeddingType = 'largeVis', M=M, perplexity=30, gamma=1/M)
r$plotEmbedding(type='PCA', show.legend=FALSE, mark.groups=TRUE, min.cluster.size=50, shuffle.colors=FALSE, font.size=3, alpha=0.3, title='clusters (largeVis)', plot.theme=theme_bw() + theme(plot.title = element_text(hjust = 0.5)))
r$getEmbedding(type='PCA', embeddingType='tSNE', perplexity=50, verbose=FALSE)
r$plotEmbedding(type='PCA', embeddingType='tSNE', show.legend=FALSE, mark.groups=TRUE, min.cluster.size=1, shuffle.colors=FALSE, font.size=3, alpha=0.3, title='clusters (tSNE)', plot.theme=theme_bw() + theme(plot.title = element_text(hjust = 0.5)))

df_allele = fread('../numbat/SRR12397038/SRR12397038_allele_counts.tsv.gz')
#out = run_numbat(count_mat = counts,
##                df_allele = df_allele,
#                 ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
#                 genome = "hg38",
#                 t = 1e-5,
#                 plot = TRUE,
#                 out_dir = '../numbat/SRR12397038/'
#)

nb = Numbat$new(out_dir = '../numbat/SRR12397038/')




plist = list()
muts = c('1a', '3b', '22b')
cnv_type = nb$joint_post %>% distinct(seg, cnv_state) %>% {setNames(.$cnv_state, .$seg)}
for (mut in muts) {
  
  plist[[mut]] = pagoda$plotEmbedding(
    alpha=0.8,
    size=1, 
    plot.na = F, 
    colors = nb$joint_post %>%
      filter(seg == mut) %>%
      {setNames(.$p_cnv, .$cell)},
    show.legend = T,
    mark.groups = F,
    plot.theme = theme_bw(),
    title = paste0(mut, '(', cnv_type[muts], ')')
  ) +
    scale_color_gradient2(low = 'royalblue', mid = 'white', high = 'red3', midpoint = 0.5, limits = c(0,1), name = 'Posterior')
}
wrap_plots(plist, guides = 'collect')