library(phyloseq)
library(vegan)
library(DESeq2)
library(ggplot2)
library(shiny)
library(lubridate)
library(tidyr)
library(dplyr)
library(viridis)
library(grid)
# library(rdrop2)
library(xtable)
library(pheatmap)
library(grid)


functionpath <- "./Functions"
source(file.path(functionpath, "_n_000_helper_functions.R"))
source(file.path(functionpath, "TeachingFunctions.R"))



# datapath <- "/Users/jvb740/MarieCurie_Work/Teaching/20180920_MicrobiomePipelineShiny_CourseLarsJuhlJensen/20180821_SinyApp_TeachingCRC/TeachingApp2/"
datapath <- "."

rv <- list()

group_var <- "CRC_status"
group_var_levels <- c("control", "case") # defines the order of the groups in all plots. If set to NULL:
color_levels <- c(cbPalette[4], cbPalette[7])
names(color_levels) <- group_var_levels

# - get directly to the pheatmap problem -
rv$ps <- readRDS(file = file.path(datapath, "phyloseq_object_Zeller_CRC.rds"))
min_obs <- 0L
prevalence <- 50
rv$ps_filt <- phyloseq::filter_taxa(rv$ps, function(x){(sum(x > min_obs) > (prevalence/100)*length(x))}, prune = TRUE)
plot_ps <- rv$ps_filt
zero_color <- "gray"
sample_colors <- list(color_levels)
names(sample_colors) <- group_var

max_abundance_for_color <- 2000


plot_heatmap_physeq(physeq = plot_ps, sample_colors = sample_colors, taxa_info_df = NULL, taxa_colors = NULL, 
                    taxa_annotation = NULL, max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                    zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = FALSE, drop_color_levels = TRUE,
                    border_color = NA, 
                    cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                    annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 16, 
                    fontsize_row = 8, fontsize_col = 8, fontsize_number = 12)

assign_variables_from_function_call(list(physeq = plot_ps, sample_colors = sample_colors, taxa_info_df = NULL, taxa_colors = NULL, 
                                     taxa_annotation = NULL, max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                                     zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = FALSE, drop_color_levels = TRUE,
                                     border_color = NA, 
                                     cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                                     annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 16, 
                                     fontsize_row = 8, fontsize_col = 8, fontsize_number = 12))


# --




# - Tab2: explore -
rv$ps <- readRDS(file = file.path(datapath, "phyloseq_object_Zeller_CRC.rds"))

OTU <- t(as(otu_table(rv$ps), "matrix"))
OTU <- round(OTU)


psShow <- capture.output(rv$ps)
psShow <- as.data.frame(psShow, nrow = 4)
psShow <- psShow[2:4, , drop = FALSE]
colnames(psShow) <- "object component"
psShow <- tidyr::separate(psShow, col = "object component", into = c("object component", "dimension"), sep = ":")

# -- calcPhyla --
Phyla <- check_phyla_distribution_NA(physeq = rv$ps)
PhylaForColor <- check_phyla_distribution(rv$ps)
if (nrow(PhylaForColor) < 16) {
        phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), QuantColors15)
} else {
        phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), c(QuantColors15, viridis(nrow(PhylaForColor)-15)))
}
rv$phylum_colors <- phylum_colors
# ----

# -- plot samples --
group_var <- "CRC_status"
group_var_levels <- c("control", "case") # defines the order of the groups in all plots. If set to NULL:
color_levels <- c(cbPalette[4], cbPalette[7])
names(color_levels) <- group_var_levels

Tr <- plot_sample_bars(physeq = rv$ps, x = "Sample", y = "Abundance", group_var, color_levels, fill = "Phylum",
                             color_sample_names = TRUE, col_vec = rv$phylum_colors, facet_grid = NULL, order_by_firmicutes = FALSE)
# ----
# --
##



########
# - Tab3: library size adjustment -
# -- caluclation of Size Factors using poscounts --
rv$ps <- readRDS(file = file.path(datapath, "phyloseq_object_Zeller_CRC.rds"))

prevalence_for_sf = 70
min_obs <- 0L

ps_sf_filt <- phyloseq::filter_taxa(rv$ps, function(x){(sum(x > min_obs) > (prevalence_for_sf/100)*length(x))}, prune = TRUE)


SFs <- calc_SFs(physeq = ps_sf_filt)

SFs_poscounts <- calc_SFs_DESeq(rv$ps, type = "poscounts", group_var = group_var)

SFs_RA <- sample_sums(rv$ps)/gm_own(sample_sums(rv$ps), zeros.count = FALSE)


library_size_adjust_list <- simply_adjust_LS(rv$ps, SFs = SFs) 
ps_tca <- library_size_adjust_list[[1]]

library_size_adjust_list <- simply_adjust_LS(rv$ps, SFs = SFs_poscounts) 
ps_tca_poscounts <- library_size_adjust_list[[1]]

library_size_adjust_list <- simply_adjust_LS(rv$ps, SFs = SFs_RA) 
ps_tca_RA <- library_size_adjust_list[[1]]
# ----
# -- plotting the data --

psP <- phyloseq::tax_glom(rv$ps, taxrank = "Phylum", NArm = FALSE)
ps_tcaP <- phyloseq::tax_glom(ps_tca, taxrank = "Phylum", NArm = FALSE)
bar_plot_ps_vs_ps_tca <- plot_sample_bars_compare(physeq = psP, physeq2 = ps_tcaP, x = "Sample", y = "Abundance", group_var = group_var, color_levels = color_levels, color_sample_names = TRUE, fill = "Phylum", col_vec = phylum_colors, order_by_raw_counts = TRUE)
# ----
# --
##





########
# - Tab4: sparsity plot -
# -- caluclation of heatmap --
sample_colors <- list(color_levels)
names(sample_colors) <- group_var

max_abundance_for_color <- 2000

plot_ps <- ps_tca

p <- plot_heatmap_physeq(plot_ps, sample_colors = sample_colors, taxa_info_df = NULL, taxa_colors = NULL, 
                         taxa_annotation = NULL, max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                         zero_color = "gray", color_function = viridis, color_steps_bw_markers = 10, log_transform = FALSE, drop_color_levels = TRUE,
                         border_color = NA, 
                         cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                         annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 14, 
                         fontsize_row = 12, fontsize_col = 12, fontsize_number = 12)

grid::grid.newpage()
grid::grid.draw(p$gtable)


# ----




########
# - Tab5: filtering -
# -- testFiltering --
prevalence <- 20
TrList <- plot_ab_pev_distributions(rv$ps, prevalence = prevalence)

ps_filt <- phyloseq::filter_taxa(rv$ps, function(x){(sum(x > min_obs) > (prevalence/100)*length(x))}, prune = TRUE)

filterList <- visualize_filtering(physeq = rv$ps, prevalence = prevalence, taxa_sums_quantile = 100, phylum_colors = phylum_colors)

# ----


########
# - Tab6: beta diversity -
# -- testFiltering --
measure <- "jsd"

dist_list <- calc_beta_div_distances(ps_filt, dist_methods = measure, group_var = group_var, compare = group_var_levels)

group_factor <- sample_data(ps_filt)[[group_var]]
# make sure group_factor fits to dist_list, i.e. only keep samples covered by compare = group_var_levels!
group_factor <- factor(group_factor[group_factor %in% group_var_levels], levels = group_var_levels, ordered = T)
adonis_list <- lapply(dist_list, function(dist_obj){
        loop_vegan_adonis(dist_obj = dist_obj, group_fac = group_factor)
})


pcoas <- calc_ordination_from_distances(ps_filt, group_var = group_var, dist_list = dist_list, color_levels = color_levels, ordination_type = "PCoA", shape = NULL, coord_cor = TRUE, phylum_colors = phylum_colors) 

TrList_samples <- pcoas[["ordination_Tr_samples"]]
# ----



########
# - Tab7: fisher -
# -- testfisher --
physeq_to_test <- ps_filt
physeq_to_test <- phyloseq::transform_sample_counts(physeq_to_test, function(x){x/sum(x)})
diff_ab_df <- test_diffs_in_prevalence_single(physeq = physeq_to_test, group_var = group_var, compare = group_var_levels, p.adj.method = "fdr", minCount = 0L, symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = NULL)
taxa_hit_df <- hit_list[["hit_table"]]

significance_colors <- brewer.pal(4, "Reds")
significance_colors <- c(rev(significance_colors), "gray", "violet")
names(significance_colors) = c("****", "***", "**", "*", "ns", "?")
taxa_colors <- list("signi_adj" = significance_colors, "Phylum" = phylum_colors)
sample_colors <- list(color_levels)
names(sample_colors) <- group_var
taxa_annotation <- taxa_hit_df$Annotation

p <- plot_heatmap_physeq(physeq_to_test, sample_colors = sample_colors, taxa_info_df = head(taxa_hit_df, 40), taxa_colors = taxa_colors, 
                         taxa_annotation = head(taxa_annotation, 40), max_abundance_for_color = .08, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                         zero_color = "gray", color_function = viridis, color_steps_bw_markers = 10, log_transform = FALSE, drop_color_levels = TRUE,
                         border_color = NA, 
                         cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                         annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 14, 
                         fontsize_row = 12, fontsize_col = 12, fontsize_number = 12)

# --

# - get a more informative taxa annotation of the hits -
taxa_annotation <- taxa_hit_df$Annotation

# ----


tt <- as.data.frame(tax_table(ps_phylum))[,2:7]

### Phylum part 
# - -
ps_phylum <- phyloseq::tax_glom(ps_filt, taxrank = "Phylum", NArm = FALSE)

taxonomic_level <- "Phylum"

# --

taxa_annotation <- get_taxon_names(as.data.frame(tax_table(ps_phylum))[,2:ncol(tax_table(ps_phylum))])
taxa_annotation <- make.unique(taxa_annotation)
taxa_order <- c(names(phylum_colors), taxa_annotation[!taxa_annotation %in% names(phylum_colors)])


# - plot Firmicutes to all other phyla ratio plots NB: you could change taxa_den to maybe only Bacteroidetes -
# NB: if you would like the order based on pValues/significance choose tax_order = NULL
FirmicutesRatioPlots <- plot_taxa_ratios_AllLevels(physeq = ps_phylum, group_var = group_var, color_levels = color_levels, tax_names = taxa_annotation, taxa_nom = "Firmicutes", taxa_den = NULL, tax_order = taxa_order, test = "wilcox.test", p_adjust_method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

# --





# - calculate ratio matrixes  -
# NB: calculates Matrixes including all samples in physeq (alternative use calculate_raw_TbTmatrixesSingle)
ps_phylum <- phyloseq::tax_glom(ps_filt, taxrank = "Phylum", NArm = FALSE)

taxonomic_level <- "Phylum"
raw_TbTmatrixes <- calculate_raw_TbTmatrixes(physeq = ps_phylum)
# raw_TbTmatrixes <- lapply(raw_TbTmatrixes, log10)
taxa_annotation <- get_taxon_names(as.data.frame(tax_table(ps_phylum))[,2:7])
taxa_annotation <- make.unique(taxa_annotation)
taxa_order <- c(names(phylum_colors), taxa_annotation[!taxa_annotation %in% names(phylum_colors)])
# --


TbT_tile <- create_raw_TbT_TilePlot(TbTmatrixes = raw_TbTmatrixes, physeq = ps_phylum, group_var = group_var, color_levels = color_levels, signi_level = 0.05, tax_names = taxa_annotation, tax_order = taxa_order, test = "wilcoxon", p_adjust_method = "none")





