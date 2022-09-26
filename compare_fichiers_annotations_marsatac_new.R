old = readRDS("/home/rparmentier/Bureau/sc_ATAC_analysis/Git_sc_ATAC_analysis/data/Annotation_TSS_pm1kb_int_ex_53utr_ctcf_cpg_woThisto_FANTOM5_prom_gr.rds")
old_df = as.data.frame(old) 
old_sf_prom = old_df %>% dplyr::filter(annotation == "FANTOM5_promoter") %>% dplyr::select(seqnames, start, end)
old_cpg = old_df %>% dplyr::filter(annotation == "CpG Island") %>% dplyr::select(seqnames, start, end)


new_df = as.data.frame(hg19_annots_gr)
df_df_fantom = new_df %>% dplyr::filter(type == "hg19_enhancers_fantom") %>% dplyr::select(seqnames, start, end)
new_df_prom = new_df %>% dplyr::filter(type == "hg19_genes_promoters") %>% dplyr::select(seqnames,start,end)
new_df_cpg = new_df %>% dplyr::filter(type == "hg19_cpg_islands") %>% dplyr::select(seqnames,start,end)

# Rien en commun ! d'où sort le fichier qu'on utilisait pour les annotations ???
a = inner_join(old_sf_prom, df_df_fantom)
b = inner_join(old_sf_prom, new_df_prom)
cpg = inner_join(new_df_cpg, old_cpg)


http://bioconductor.org/packages/devel/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub.html
library(AnnotationHub)
ah = AnnotationHub()
dm = query(ah, c("Grange", "Homo sapiens"))
dm
df = mcols(dm)
df
