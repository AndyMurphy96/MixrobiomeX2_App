# --
#######################################
### FUNCTION: check_phyla_distribution#
#######################################
# NB: throws error if "Phylum" is not in colnames(tax_table(physeq))
# outputs a data.frame, summarising the taxa distribution over the different phyla.
## INPUT:
# physeq: physeq object
## OUTPUT:
# data.frame, summarising the taxa distribution over the different phyla.
# the Phyla are ordered so the Phylum with most counts (PC_of_counts) is on top, no of taxa is used to break ties in the ordering
# columns: PC stands for percentage, 
# mean/median_taxa_sum are the mean/median of the taxa_sums (total counts over all samples) of the taxa in the respective phylum
# other columns should be clear

check_phyla_distribution <- function(physeq) {
        
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0))
        
        
        df_ab_prev <- cbind(df_ab_prev, as.data.frame(unclass(tax_table(physeq))))
        
        PhylaDistribution <- dplyr::summarise(group_by(df_ab_prev, Phylum), 
                                              taxa = n(), 
                                              PC_of_taxa = round(100*taxa/ntaxa(physeq),1),
                                              PC_of_counts = round(100*sum(total_counts)/sum(otu_table(physeq)), 1),
                                              PC_of_prevalence = round(100*sum(prevalence)/sum(otu_table(physeq) != 0), 1),
                                              mean_taxa_sum = round(mean(total_counts)),
                                              median_taxa_sum = round(median(total_counts)),
                                              mean_prevalence_in_PC = round(100*mean(prevalence)/nsamples(physeq), 1)) %>% 
                arrange(desc(PC_of_counts), desc(taxa), desc(PC_of_prevalence))
        
        PhylaDistribution
        
}
# --



# --
#######################################
### FUNCTION: check_phyla_distribution_NA#
#######################################
# NB: throws error if "Phylum" is not in colnames(tax_table(physeq))
# outputs a data.frame, summarising the taxa distribution over the different phyla.
## INPUT:
# physeq: physeq object
## OUTPUT:
# data.frame, summarising the taxa distribution over the different phyla.
# the Phyla are ordered so the Phylum with most counts (PC_of_counts) is on top, no of taxa is used to break ties in the ordering
# columns: PC stands for percentage, 
# mean/median_taxa_sum are the mean/median of the taxa_sums (total counts over all samples) of the taxa in the respective phylum
# other columns should be clear

check_phyla_distribution_NA <- function(physeq) {
        
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0))
        
        
        df_ab_prev <- cbind(df_ab_prev, as.data.frame(unclass(tax_table(physeq))))
        df_ab_prev$Phylum <- as.character(df_ab_prev$Phylum)
        df_ab_prev$Phylum[is.na(df_ab_prev$Phylum)] <- paste(df_ab_prev$Kingdom[is.na(df_ab_prev$Phylum)], "_NA", sep = "")
        
        PhylaDistribution <- dplyr::summarise(group_by(df_ab_prev, Phylum), 
                                              taxa = n(), 
                                              PC_of_taxa = round(100*taxa/ntaxa(physeq),1),
                                              PC_of_counts = round(100*sum(total_counts)/sum(otu_table(physeq)), 1),
                                              PC_of_prevalence = round(100*sum(prevalence)/sum(otu_table(physeq) != 0), 1),
                                              mean_taxa_sum = round(mean(total_counts)),
                                              median_taxa_sum = round(median(total_counts)),
                                              mean_prevalence_in_PC = round(100*mean(prevalence)/nsamples(physeq), 1)) %>% 
                arrange(desc(PC_of_counts), desc(taxa), desc(PC_of_prevalence))
        
        PhylaDistribution
        
}
# --





# --
#######################################
#### plot_sample_bars
#######################################
# based on plot_bar from phyloseq, difference, orders and colors the Samples based on group_var, and orders the fill based on abundance
# I guess inputs can be guessed on

plot_sample_bars <- function(physeq, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                             color_sample_names = TRUE, col_vec = NULL, facet_grid = NULL, order_by_firmicutes = TRUE){
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the phyloseq object.")
        }
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all names in names(color_levels)are actually levels in the group_var column.")
        }
        
        
        # in case you do not want to see all samples
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
        }
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
        
        # if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- phyloseq::psmelt(physeq)
        
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA" # NB: pools all AN
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        
        # order samples according to levels or Firmicutes
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        
        if (order_by_firmicutes) {
                mdf_firmicutes <- dplyr::filter(mdf, Phylum == "Firmicutes") %>% arrange_(group_var, "Abundance")
                mdf$Sample <- factor(mdf$Sample, levels = mdf_firmicutes$Sample, ordered = TRUE)
        } else {
                mdf$Sample <- factor(mdf$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        }
        
        
        
        # - define names of x axis using color_levels (which must be a named character vector) -        
        colxaxis <- color_levels[LookUpDF$Group]
        # --
        
        if (is.null(col_vec)){
                if (length(levels(mdf[, fill])) <= 15) {
                        fill_colors <- make_color_vector(mdf[, fill], rev(QuantColors15[1:length(levels(mdf[, fill]))]))
                } else {
                        fill_colors <- make_color_vector(mdf[, fill], viridis(length(levels(mdf[, fill]))))
                }
                
        } else {
                fill_colors <- col_vec
        }
        
        
        
        Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
        Tr <- Tr + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors) +
                xlab("") +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        if (!is.null(facet_grid)) {
                formulation <- as.formula(paste("~ ", facet_grid, sep = ""))
                Tr <- Tr + facet_grid(formulation)
        }
        
        
        Tr
}
# --





# --
#######################################
#### plot_sample_bars_compare
#######################################
# generates abundance barplots (see plot_bar_own) to compare ps to ps_tca, i.e. to see how SFs adjustment affects the abundances 

plot_sample_bars_compare <- function(physeq, physeq2, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                                     color_sample_names = TRUE, col_vec = NULL, order_by_raw_counts = TRUE){
        
        # - prepare mdf of ps physeq -
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        # in case you do not want to see all samples
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
        }
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
        
        # if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- phyloseq::psmelt(physeq)
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA" # NB: pools all NA
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        
        # --
        
        # - prepare mdf of ps_tca physeq -
        if(taxa_are_rows(physeq2)) { physeq2 <- t(physeq2) }
        
        if (!all(unique(sample_data(physeq2)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq2)[sample_data(physeq2)[[group_var]] %in% names(color_levels)]
                physeq2 <- prune_samples(samples = keepSamples, physeq2)
                sample_data(physeq2)[[group_var]] <- factor(sample_data(physeq2)[[group_var]], levels = names(color_levels), order = TRUE)
        }
        
        if (!is.factor(sample_data(physeq2)[[group_var]])) {sample_data(physeq2)[[group_var]] <- factor(sample_data(physeq2)[[group_var]], levels = names(color_levels), order = TRUE)}
        
        
        mdf2 <- phyloseq::psmelt(physeq2)
        
        # order fill levels according to abundance over all samples
        mdf2[, fill] <- as.character(mdf2[, fill])
        mdf2[is.na(mdf2[, fill]), fill] <- "NA"
        # sums <- group_by_(mdf2, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf2[, fill] <- factor(mdf2[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        # --
        
        mdf$Typer <- "before"
        mdf2$Typer <- "after SF adjustment"
        
        mdf <- rbind(mdf, mdf2)
        mdf$Typer <- factor(mdf$Typer, levels = c("before", "after SF adjustment"), ordered = TRUE)
        
        
        # order samples according to group_var levels and potentially by total counts
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        
        if (order_by_raw_counts) {
                rawSampleSizes <- data.frame("Sample" = sample_names(physeq), "Group" = sample_data(physeq)[[group_var]], "Total" = sample_sums(physeq))
                rawSampleSizes <- arrange(rawSampleSizes, Group, Total)
                mdf$Sample <- factor(mdf$Sample, levels = rawSampleSizes$Sample, ordered = TRUE)
                
        } else {
                mdf$Sample <- factor(mdf$Sample, levels = LookUpDF$Sample, ordered = TRUE)
                
        }
        
        
        # - define names of x axis using color_levels (which must be a named character vector) -        
        colxaxis <- color_levels[LookUpDF$Group]
        # --
        
        
        if (is.null(col_vec)){
                if (length(levels(mdf[, fill])) <= 15) {
                        fill_colors <- make_color_vector(mdf[, fill], rev(QuantColors15[1:length(levels(mdf[, fill]))]))
                } else {
                        fill_colors <- make_color_vector(mdf[, fill], viridis(length(levels(mdf[, fill]))))
                }
                
        } else {
                fill_colors <- col_vec
        }
        
        
        Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
        Tr <- Tr + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors, na.value = "red") +
                xlab("") +
                facet_wrap(~ Typer, ncol = 1, scales = "free_y") +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        Tr
}
# --




# --
#######################################
### calc_SFs
#######################################
# implementation of DESeq2 library size similar to estimateSizeFactorsForMatrix
# in difference to adjust_LS (obsolete) ignore.zero.ratios is always TRUE (so 0 ratios are always ignored)


calc_SFs <- function(physeq, zeros.count = FALSE, percentile = 50)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (taxa_are_rows(physeq)) {
                SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } else {
                SFs <- apply(as(otu_table(physeq), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } 
        
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        SFs
}
# --



# --
#######################################
### calc_SFs_DESeq
#######################################


calc_SFs_DESeq <- function(physeq, group_var, type = "poscounts"){
        
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], ordered = FALSE)
        
        DES = phyloseq::phyloseq_to_deseq2(physeq, formula(paste("~", group_var)))
        
        dds <- estimateSizeFactors(DES, type = type)
        
        sizeFactors(dds)
        
        
}

# --




# --
#######################################
### simply_adjust_LS
#######################################
# implementation of DESeq2 library size similar to estimateSizeFactorsForMatrix
# in difference to adjust_LS (obsolete) ignore.zero.ratios is always TRUE (so 0 ratios are always ignored)

## Input:
# physeq
# zeros.count: zeros.count of gm_own, if TRUE 0 will be considered when calculating the geometric mean
# if FALSE not and thus the geometric means will be bigger (see gm_own)
# percentile: the percentile to select the size factor SF of a sample based on its count ratios to the reference sample. In
# DESeq percentile = 50, i.e. stats::median is used. 
# plots: if TRUE SFs and plots will be given in an output list, otherwise just the adjusted physeq is returned


simply_adjust_LS <- function(physeq, SFs = NULL, zeros.count = FALSE, percentile = 50)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors  unless given --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (is.null(SFs)){
                
                if (taxa_are_rows(physeq)) {
                        SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                        SFs <- SFs/exp(mean(log(SFs)))
                } else {
                        SFs <- apply(as(otu_table(physeq), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                        SFs <- SFs/exp(mean(log(SFs)))
                }
                
        }
        
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        
        # --- 3: calculate the new counts and put into a physeq object
        
        if(taxa_are_rows(physeq)){
                if (!identical(colnames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),2,SFs, "/")
                phynew <- physeq
                otu_table(phynew) <- otu_table(Mat, taxa_are_rows = TRUE)
        } else {
                if (!identical(rownames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),1,SFs, "/") 
                phynew <- physeq
                otu_table(phynew) <- otu_table(Mat, taxa_are_rows = FALSE)
        }
        
        
        
        list(physeq = phynew, SFs = SFs)
        
}
# --



# --
#########################
## plot_heatmap_physeq ##
#########################
# function that uses pheatmap to draw a heatmap of the physeq object, or a pruned version of it determined by taxa_info_df
# INPUTs
# physeq = physeq object
# sample_colors = named list of named color character vectors. The sample_data(physeq) data.frame is used as sample_info_df to color the samples (columns) in the heatmap.
# Only those variables/columns in sample_info_df that are names in sample_colors are considered for coloring the samples. I.e. if sample_colors = NULL no coloring occurs.
# NB: if a variable is a name in sample_colors but the color_character is NA, default colors will be used.
# NB!!: The named character vectors in sample_colors are also used to prune the samples!! This works as follows:
# The character vectors in sample_colors are checked in order. If the first character vector (sample_colors[[1]]) is for example "Group" and it contains only colors for "A" and "B", 
# but sample_info_df[["Group"]] contains also "C", the C samples will be removed! The procedure continues for all entries in sample_colors
# finally the order of names(sample_colors) defines the order of the color bars, the first entry being drawn closest to the heatmap and so on. 
# taxa_info_df = a data frame similar to sample_info_df (see sample_colors) that together with taxa_colors will be used to label the taxa.
# NB: MORE importantly, it is also used to restrict the included taxa and to order the taxa!! The rownames(taxa_info_df) must be the taxa_names in physeq 
# of the taxa you want to include in the heatmap in the order you want the taxa to be shown! IF NULL all taxa will be plotted in the order given by physeq
# taxa_colors: see sample_colors
# taxa_annotation: if not NULL will be used as labels_row in pheatmap, i.e. a different labelling of the taxa instead of taxa_names(physeq) or rownames(taxa_info_df)
# max_abundance_for_color: the maximum abundance value to which the maximum color is assigned, all abundances above get also this color, if NULL max(otu_table(physeq))
# is used
# gradient_steps:  gradient_steps: vector of numbers between 0 and 1, default c(0.15, 0.3, 0.45, 1). Defines how the main colors in the gradient will be positioned between min_in_data and
# max_abundance_for_color. min_in_data (= lowest non-zero count/abundance) will be added to first position after gradient_steps have been normalized with max_abundance_for_color
# zero_color: defines the color of the 0 values (minimum Value in otu_table(physeq).
# color_function: a function such as the default viridis function that generates color strings when given a number.
# color_steps_bw_markers: numeric value defining how many equally distributed breaks will be introduced between the markers defined by gradient_steps.
# log_transform: option to log10 transform the counts
# drop_color_levels: if TRUE, color levels in sample_colors and taxa_colors are removed (not shown in legend so) when not present in the data
# rest pheatmap arguments

plot_heatmap_physeq <- function (physeq, sample_colors = NULL, taxa_info_df = NULL, taxa_colors = NULL, taxa_annotation = NULL, 
                                 max_abundance_for_color = NULL, gradient_steps = c(0.15, 0.3, 0.45, 1), zero_color = "white",
                                 color_function = viridis,
                                 color_steps_bw_markers = 10, log_transform = FALSE, drop_color_levels = FALSE,
                                 border_color = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = TRUE, annotation_names_col = TRUE, annotation_legend = TRUE,
                                 legend = TRUE, font_size = 10, fontsize_row = 8, fontsize_col = 8, fontsize_number = 6, filename = NA, ...) {
        
        # - keep only taxa in taxa_info_df in physeq -
        if (!is.null(taxa_info_df)) {
                
                if (!all(rownames(taxa_info_df) %in% taxa_names(physeq))) {
                        stop("not all taxa in taxa_info_df are taxa in physeq!")
                }
                pruned_ps <- prune_taxa(rownames(taxa_info_df), physeq)
        } else {
                pruned_ps <- physeq
                taxa_colors <- NULL
        }
        # --
        
        # - test sample_colors input and prepare sample_info_df -
        if (!is.null(sample_colors)) {
                
                if (!is.list(sample_colors)) {
                        stop("sample_colors must be a named list.")
                }
                
                
                if (!all(names(sample_colors) %in% colnames(sample_data(pruned_ps)))) {
                        stop("all names of the sample_colors list must be variables/columns in sample_data(physeq).")
                }
                
                sample_info_df <- as(sample_data(pruned_ps), "data.frame")
                sample_info_df <- sample_info_df[, match(names(sample_colors), colnames(sample_info_df)), drop = FALSE] # the match construct makes sure that sample_info_df columns are ordered after names(sample_colors)
                # which makes sure that the column coloring occurs in the whished for order
                
                # -- test all entries in sample_colors, prune pruned_ps and order sample_info_df, and set default color where necessary --
                for (i in 1:length(sample_colors)){
                        variable_name <- names(sample_colors)[i]
                        color_character <- sample_colors[[i]]
                        
                        if (!all(is.na(color_character))){
                                
                                # to not remove all samples from pruned_ps:
                                if (!(sum(names(color_character) %in% unique(sample_info_df[[variable_name]])) > 0)) {
                                        stop(paste("Not a single name of a color for ", variable_name, " matched to an entry in the corresponding column in sample_data(physeq)", sep = ""))
                                }
                                
                                # remove all samples that no color was assigned to in sample_colors, i.e.option to restrict the comparison to defined samples
                                if (!all(unique(sample_info_df[[variable_name]]) %in% names(color_character))) {
                                        
                                        keepSamples <- sample_names(pruned_ps)[sample_info_df[[variable_name]] %in% names(color_character)]
                                        pruned_ps <- prune_samples(keepSamples, pruned_ps)
                                        sample_info_df <- as(sample_data(pruned_ps), "data.frame")
                                        sample_info_df <- sample_info_df[, colnames(sample_info_df) %in% names(sample_colors), drop = FALSE]
                                }
                                
                                # remove unused color levels if asked for it
                                if (drop_color_levels && !all(names(color_character) %in% unique(sample_info_df[[variable_name]]))) {
                                        color_character <- color_character[names(color_character) %in% unique(sample_info_df[[variable_name]])]
                                        sample_colors[[i]] <- color_character
                                }
                                
                                if (!all(areColors(color_character))) {
                                        warning(paste("Not all entries in sample_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(sample_info_df, variable_name)
                                        sample_colors[[i]] <- color_character
                                }
                                
                                # define the order of the samples based on the sample_color vectors (actual ordering below with dplyr::arrange)
                                sample_info_df[[variable_name]] <- factor(sample_info_df[[variable_name]], levels = names(color_character), ordered = TRUE)
                                
                        } else {
                                color_character <- assign_default_colors(sample_info_df, variable_name)
                                sample_colors[[i]] <- color_character
                        }
                }
                sample_info_df$IDSaver <- rownames(sample_info_df)
                sample_info_df <- dplyr::arrange_(sample_info_df, names(sample_colors))
                rownames(sample_info_df) <- sample_info_df$IDSaver
                sample_info_df <- sample_info_df[, -ncol(sample_info_df), drop = FALSE] # to remove IDSaver
                # ----
        } else {
                sample_info_df <- NULL
        }
        # --
        
        # - set the taxa colors (and pool sample_colors and taxa_colors into annotation_colors) -
        if (!is.null(taxa_colors)) { # see above, guarantees also that taxa_info_id is not NULL
                
                if (!is.list(taxa_colors)) {
                        stop("taxa_colors must be a named list.")
                }
                
                if (!all(names(taxa_colors) %in% colnames(taxa_info_df))) {
                        stop("all names of the taxa_colors list must be variables/columns in taxa_info_df.")
                }
                
                taxa_info_df <- taxa_info_df[, match(names(taxa_colors), colnames(taxa_info_df)), drop = FALSE] # again match to get coloring ordered after taxa_colors
                
                # -- test all entries in taxa_colors --
                for (i in 1:length(taxa_colors)){
                        variable_name <- names(taxa_colors)[i]
                        color_character <- taxa_colors[[i]]
                        
                        if (!all(is.na(color_character))){
                                
                                
                                if (!all(unique(taxa_info_df[[variable_name]]) %in% names(color_character))) {
                                        warning(paste("There were levels in taxa_info_df at ", variable_name, " for which no color was assigned in taxa_colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(taxa_info_df, variable_name)
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                # remove unused color levels if asked for it
                                if (drop_color_levels && !all(names(color_character) %in% unique(taxa_info_df[[variable_name]]))) {
                                        color_character <- color_character[names(color_character) %in% unique(taxa_info_df[[variable_name]])]
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                
                                if (!all(areColors(color_character))) {
                                        warning(paste("Not all entries in taxa_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(taxa_info_df, variable_name)
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                # define the order of the taxa coloring based on the taxa_color vectors
                                # taxa_info_df[[variable_name]] <- factor(taxa_info_df[[variable_name]], levels = names(color_character), ordered = TRUE) # was unnecessary, pheatmap orders based on annotation_colors
                                
                        } else {
                                color_character <- assign_default_colors(taxa_info_df, variable_name)
                                taxa_colors[[i]] <- color_character
                        }
                }
                # ----
        } 
        
        annotation_colors <- c(sample_colors, taxa_colors)
        # -- 
        
        # - check or set taxa_annotation -
        if (is.null(taxa_annotation)){
                taxa_annotation <- taxa_names(pruned_ps) # or you could use 
                # taxa_annotation <- get_taxon_names(as.data.frame(tax_table(pruned_ps)))
        }
        
        if (length(taxa_annotation) != ntaxa(pruned_ps)) {
                warning("taxa_annotation did not fit in length to nrow(taxa_info_df) or ntaxa(physeq), used taxa_names")
                taxa_annotation <- taxa_names(pruned_ps)
        }
        
        taxa_annotation <- as.character(taxa_annotation)
        
        taxa_annotation[is.na(taxa_annotation)] <- "NA"
        
        taxa_annotation <- make.unique(taxa_annotation)
        # --
        
        # - generate count data frame in which taxa are rows in the order determined by taxa_info_df -
        if (taxa_are_rows(pruned_ps)) {
                pruned_ps <- t(pruned_ps)
        }
        
        DF_CT <- as(otu_table(pruned_ps), "matrix")
        DF_CT <- as.data.frame(t(DF_CT))
        
        if (!is.null(taxa_info_df)){
                DF_CT <- DF_CT[rownames(taxa_info_df), ]
        }
        # --
        
        # - order the samples in DF_CT based on sample_colors list (see above for sample_info_df) -
        if (!is.null(sample_info_df)){
                DF_CT <- DF_CT[, rownames(sample_info_df)]
        }
        # --
        
        # - test and adjust max_abundance_for_color -
        if (is.null(max_abundance_for_color)) {
                max_abundance_for_color <- max(DF_CT)
        }
        if (max_abundance_for_color < min(DF_CT) && max_abundance_for_color > max(DF_CT)) {
                max_abundance_for_color <- max(DF_CT)
        }
        # --
        
        # - set breaks and colors for pheatmap and do a possible log transform -
        ZeroValue <- min(DF_CT) # should be 0 in most cases of microbiome data, unless pseudocount had been added. In rare cases of high taxonomy it might be a non-zero count so actually min_in_data 
        min_in_data <- min(DF_CT[DF_CT > ZeroValue]) # the lowest non-zero value
        max_in_data <- max(DF_CT)
        
        # -- make sure the last entry in gradient steps = 1 --
        if (!all(gradient_steps >= 0 & gradient_steps <= 1)) {
                gradient_steps <- c(0.15, 0.3, 0.45, 1)
        }
        
        if (gradient_steps[length(gradient_steps)] != 1) {
                gradient_steps <- c(gradient_steps, 1)
        }
        # ----
        
        # you want that the final color gradient covers the values min to max_abundance_for_color
        # 0 values will be set to a different color (usually red or white), values above max_abundance_for_color should be all max_color
        # normalise gradient steps with max_abundance_for_color
        myBreaks <- gradient_steps * max_abundance_for_color
        myBreaks <- myBreaks[myBreaks > min_in_data]
        # add break at min_in_data
        myBreaks <- c(min_in_data, myBreaks)
        
        # now myBreaks goes from min_in_data up to max_abundance_for_color (provided the last gradient_steps was 1)
        myColors = color_function(length(myBreaks)) # see help pheatmap, breaks should be 1 element longer than color, now it is same legnth
        
        # myColors contains now the colors that represent the gradient_steps values (= markers).
        # now we want to introduce breaks between these markers and make linear color gradients between the markers. The number of breaks between the markers is defined
        # by color_steps_bw_markers
        
        myColors <- unlist(lapply(1:(length(myColors)-1), function(i) {
                colorRampPalette(colors = c(myColors[i], myColors[i+1]))(color_steps_bw_markers)
        }))
        
        
        
        # -- do a possible log transform using min_in_data/5 as pseudocounts --
        if (log_transform) {
                if (ZeroValue > 0){ # > 0 if pseudocount had been added or simply no zero counts in data, < 0 if data h
                        pseudocount <- ZeroValue
                } else if (ZeroValue < 0){
                        stop("you asked for log_transform but the lowest count in your data is already below 0!")
                } else {
                        pseudocount <- min_in_data/2
                }
                DF_CT[DF_CT == ZeroValue] <- pseudocount
                DF_CT <- log10(DF_CT)
                
                myBreaks <- log10(myBreaks)
                myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
                        seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
                })
                myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)]) #length(myBreaks) is now length(myBreaks) * color_steps_bw_markers + 1, the markers are at positions 1, 1+color_steps_bw_markers, 1+2*color_steps_bw_markers 
                
                if (log10(max_in_data) > myBreaks[length(myBreaks)]) {
                        myBreaks <- c(log10(pseudocount), myBreaks, log10(max_in_data))
                        myColors <- c(zero_color, myColors, myColors[length(myColors)])
                } else {
                        myBreaks <- c(log10(pseudocount), myBreaks)
                        myColors <- c(zero_color, myColors)
                }
                
        } else {
                myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
                        seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
                })
                myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)])
                
                # add max_in_data values to myBreaks, otherwise all values above max_abundance_for_color will be white
                if (max_in_data > myBreaks[length(myBreaks)]) {
                        myBreaks <- c(ZeroValue, myBreaks, max_in_data)
                        myColors <- c(zero_color, myColors, myColors[length(myColors)])
                } else {
                        myBreaks <- c(ZeroValue, myBreaks)
                        myColors <- c(zero_color, myColors)
                }
        }
        # ----
        # --
        
        # return()
        # I know it works until here
        
        # - not necessary, just for clarity -
        if (is.null(annotation_colors)) {annotation_colors <- NA}
        if (is.null(sample_info_df)) {sample_info_df <- NA}
        if (is.null(taxa_info_df)) {taxa_info_df <- NA}
        # --
        
        # --
        hm.parameters <- list(DF_CT, color = myColors, breaks = myBreaks, border_color = border_color, cluster_cols = cluster_cols, cluster_rows = cluster_rows,
                              show_rownames = show_rownames, show_colnames = show_colnames, annotation_col = sample_info_df,
                              annotation_row = taxa_info_df, annotation_colors = annotation_colors, labels_row = taxa_annotation, annotation_names_row = annotation_names_row,
                              annotation_names_col = annotation_names_col, annotation_legend = annotation_legend, legend = legend, font_size = font_size,
                              fontsize_row = fontsize_row, fontsize_col = fontsize_col, fontsize_number = fontsize_number)
        do.call("pheatmap", hm.parameters)
        # --
        
        
}
# --





# --
#######################################
#### plot_ab_pev_distributions
#######################################


plot_ab_pev_distributions <- function(physeq, prevalence = 5) {
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        seqtab <- as(otu_table(physeq), "matrix")
        
        FinalNumbersSeq <- data.frame(Sequence = colnames(seqtab), InNumberSamples = colSums(seqtab != 0), TotalCounts = colSums(seqtab))
        # a look at it shows you that seqtab is ordered after total counts
        FinalNumbersSeq <- group_by(FinalNumbersSeq, InNumberSamples)
        CountDistribution <- dplyr::summarise(FinalNumbersSeq, No_ASVs = n(), TotalCounts = sum(TotalCounts))
        CountDistribution$CumSumUnique <- rev(cumsum(rev(CountDistribution$No_ASVs)))
        CountDistribution$CumPerCUnique <- rev(cumsum(rev(CountDistribution$No_ASVs/ncol(seqtab))))
        CountDistribution$CumSumTotal <- rev(cumsum(rev(CountDistribution$TotalCounts)))
        CountDistribution$CumPerCTotal <- rev(cumsum(rev(CountDistribution$TotalCounts/sum(colSums(seqtab)))))
        
        PCValue <- ceiling((prevalence/100)*dim(seqtab)[1]) # tells you in how many samples a SV must be present to meet the prevalence 
        
        # Diff <- CountDistribution$InNumberSamples - PCValue
        # index <- which.max(Diff[Diff<0]) + which.min(Diff[Diff>=0])
        index <- which(CountDistribution$InNumberSamples >= PCValue)[1]
        PCKeptAtPCValue <- CountDistribution$CumPerCTotal[index]
        SVskeptAtPCValue <- CountDistribution$CumPerCUnique[index]
        
        # The number of samples the SVs are present in
        Tr <- ggplot(CountDistribution, aes(x = InNumberSamples, y = No_ASVs))
        Tr <- Tr + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("number of taxa") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        
        In1Index <- which(CountDistribution$InNumberSamples == 1)
        if (length(In1Index) != 0) {
                Tr <- Tr + ggtitle(paste(CountDistribution$No_ASVs[In1Index], " of ", CountDistribution$CumSumUnique[1], " taxa (", round(100*CountDistribution$No_ASVs[In1Index]/CountDistribution$CumSumUnique[1], 1), " %)", " were only found in 1 sample", sep = ""))
        } 
        
        
        # Cumulative Percentage of SVs
        Tr1 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = CumPerCUnique))
        Tr1 <- Tr1 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of taxa") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr1 <- Tr1 + 
                geom_hline(yintercept = SVskeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", CountDistribution$CumSumUnique[index], " of ", CountDistribution$CumSumUnique[1], 
                              " taxa (", round(100*CountDistribution$CumSumUnique[index]/CountDistribution$CumSumUnique[1], 1), " %) have higher prevalence", sep = ""))
        
        
        Tr2 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = TotalCounts))
        Tr2 <- Tr2 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("total counts of taxa with given prevalence") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr2 <- Tr2 + ggtitle(paste(CountDistribution$CumSumTotal[1] - CountDistribution$CumSumTotal[index], " of ",
                                   CountDistribution$CumSumTotal[1], " (", round(100*(CountDistribution$CumSumTotal[1] - CountDistribution$CumSumTotal[index])/CountDistribution$CumSumTotal[1], 2),
                                   " %) counts are from taxa present in less than ", round((prevalence/100)*dim(seqtab)[1],1), " samples.", sep = ""))
        
        Tr3 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = CumPerCTotal))
        Tr3 <- Tr3 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of counts") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr3 <- Tr3 + 
                geom_hline(yintercept = PCKeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", CountDistribution$CumSumTotal[index], " of ", CountDistribution$CumSumTotal[1], 
                              " counts (", round(100*CountDistribution$CumSumTotal[index]/CountDistribution$CumSumTotal[1], 1), " %) would remain", sep = ""))
        
        
        TrList <- list(Tr, Tr1, Tr2, Tr3, CountDistribution)
        
        TrList
        
}
# --





# --
#######################################
#### visualize_filtering
#######################################


visualize_filtering <- function(physeq, prevalence, taxa_sums_quantile, phylum_colors = NULL){ 
        
        if (taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0),
                                 sparsity = colSums(as(otu_table(physeq), "matrix") == 0),
                                 mean_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){mean(x[x > 0])}),
                                 median_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){median(x[x > 0])})
        )
        
        df_ab_prev <- cbind(df_ab_prev, as.data.frame(unclass(tax_table(physeq))))
        
        # - adjust color and order of the phyla in the following plots - 
        colori <- "Phylum"
        
        if (!is.null(phylum_colors)){
                df_ab_prev[[colori]] <- as.character(df_ab_prev[[colori]])
                df_ab_prev[[colori]][is.na(df_ab_prev[[colori]])] <- "NA" # NB: pools NAs
                if (!all(unique(df_ab_prev[[colori]]) %in% names(phylum_colors))){
                        stop("provided phylum_colors did not cover all Phyla in physeq")
                }
                df_ab_prev[[colori]] <- factor(df_ab_prev[[colori]], levels = names(phylum_colors), ordered = TRUE)
                custom_colors <- phylum_colors
        } else {
                CountOrder <- dplyr::group_by_(df_ab_prev, colori) %>% dplyr::summarise(total_count_sum = sum(total_counts)) %>% dplyr::arrange(desc(total_count_sum))
                
                CountOrder[[colori]] <- as.character(CountOrder[[colori]])
                CountOrder[[colori]][is.na(CountOrder[[colori]])] <- "NA"
                
                if (nrow(CountOrder) <= 15){
                        custom_colors <- make_color_vector(CountOrder[[colori]], QuantColors15)
                } else {
                        custom_colors <- make_color_vector(CountOrder[[colori]], viridis(nrow(CountOrder)))
                }
                
                df_ab_prev[[colori]] <- as.character(df_ab_prev[[colori]])
                df_ab_prev[[colori]][is.na(df_ab_prev[[colori]])] <- "NA"
                df_ab_prev[[colori]] <- factor(df_ab_prev[[colori]], levels = names(custom_colors), ordered = TRUE)
                
        }
        # --
        
        prev_thresh <- (prevalence/100)*nsamples(physeq)
        abund_thresh <- quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100)
        
        df_ab_prev_filt <- dplyr::filter(df_ab_prev, prevalence > prev_thresh | total_counts > abund_thresh)
        
        no_samples <- nsamples(physeq)
        shade_df <- data.frame(total_counts = 0, prevalence = 0)
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts, y = prevalence))
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                geom_point(aes_string(col = colori), size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts (taxa_sums())") +
                scale_color_manual("", values = custom_colors) +
                theme_bw() +
                ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " taxa (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                              " %) and ", round(sum(df_ab_prev_filt$total_counts)), " of ", round(sum(df_ab_prev$total_counts)), " counts (",
                              round((sum(df_ab_prev_filt$total_counts)/sum(df_ab_prev$total_counts))*100, 1), " %) remain", sep = ""))
        
        
        
        Tr_prev_vs_log10_ab_wrap <- ggplot(df_ab_prev, aes(x = total_counts, y = prevalence))
        Tr_prev_vs_log10_ab_wrap <- Tr_prev_vs_log10_ab_wrap +
                geom_point(size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts (taxa_sums())") +
                facet_wrap(~ Phylum) +
                scale_color_manual("", values = custom_colors) +
                theme_bw() +
                ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " taxa (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                              " %) and ", round(sum(df_ab_prev_filt$total_counts)), " of ", round(sum(df_ab_prev$total_counts)), " counts (",
                              round((sum(df_ab_prev_filt$total_counts)/sum(df_ab_prev$total_counts))*100, 1), " %) remain", sep = ""))
        
        
        
        Tr_prev_vs_log10_ab_wrap <- Tr_prev_vs_log10_ab +
                facet_wrap(~ Phylum) +
                theme(legend.position = "none")
        
        out <- list(Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab,
                    Tr_prev_vs_log10_ab_wrap = Tr_prev_vs_log10_ab_wrap)
        
}
# --



# --
#######################################
### FUNCTION: calc_beta_div_distances
#######################################
# see unlist(phyloseq::distanceMethodList) for all available methods


calc_beta_div_distances <- function(physeq, dist_methods = c("bray"), group_var = NULL, compare = NULL) {
        
        
        if (!is.null(group_var) || !is.null(compare)){
                
                if(! group_var %in% colnames(sample_data(physeq))) {
                        stop("The given group_var is not a variable in the sample data of the phyloseq object.")
                }
                
                if (!all(compare %in% unique(sample_data(physeq)[[group_var]]))) {
                        stop("Not all names in compare are actually levels in the group_var column.")
                }
                
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% compare]
                physeq <- prune_samples(keepSamples, physeq)
                # not sure if taxa that are not present in a single sample after the prune affect distance?:
                # tested: at least not for bray curtis, and jsd. When you find one, uncomment:
                # physeq <- phyloseq::subset_taxa(physeq, taxa_sums(ps) != 0) 
                
        }
        
        
        dist_list <- vector("list", length(dist_methods))
        names(dist_list) = dist_methods
        
        for (i in dist_methods) {
                iDist <- phyloseq::distance(physeq, method=i) # so for "bray" same as vegan::vegdist(x = as(otu_table(physeq), "matrix"), method = "bray")
                dist_list[[i]] = iDist
        }
        
        return(dist_list)
        
}
# --







# --
#######################################
### FUNCTION: loop_vegan_adonis
#######################################
# Function is very much based on pairwise.perm.manova {RVAideMemoire}
# but it also records R2 while looping through vegan::adonis, and generates a
# result data frame in which the results are shown in the order of the group_fac levels
# INPUT:
# dist_obj: dist object
# group_fac: the factor that groups the samples
# nperm: permutations in adonis
# p.adj.method: method to adjust p.values
# symnum.args: symbols linked to p-value cutpoints
# OUTPUT:
# data.frame showing p.values and R2 and adjusted p.values for the different between group comparisons


loop_vegan_adonis <- function(dist_obj, group_fac, nperm = 999, 
                              p.adj.method = "none", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        if (!("dist" %in% class(dist_obj))){
                stop("dist_obj must be of class dist")
        }
        
        if (length(dist_obj) != length(group_fac)*(length(group_fac)-1)/2) {
                stop("The number of samples in group_fac does not fit to the number of distances in dist_obj. Check and come back.")
        }
        
        group_fac <- factor(group_fac)
        
        # - get all level combinations within group_fac -
        co = combn(levels(group_fac),2)
        # -- 
        
        p_vals <- vector(mode = "numeric", length = ncol(co))
        r2s <- vector(mode = "numeric", length = ncol(co))
        F.Models <- vector(mode = "numeric", length = ncol(co))
        pairs <- vector(mode = "character", length = ncol(co))
        
        # - loop vegan::adonis for all level combinations and generate df -
        for (k in 1:ncol(co)){
                group_fac2 <- droplevels(group_fac[group_fac %in% co[,k]])
                dist_obj_mat <- as.matrix(dist_obj)
                rows <- which(group_fac %in% levels(group_fac2))
                dist_obj2 <- as.dist(dist_obj_mat[rows, rows])
                fit <- vegan::adonis(dist_obj2 ~ group_fac2, permutations = nperm)
                p_vals[k] <- fit$aov.tab[1, "Pr(>F)"]
                r2s[k] <- fit$aov.tab[1, "R2"]
                F.Models[k] <- fit$aov.tab[1, "F.Model"]
                pairs[k] <- paste(co[,k], collapse = " vs ")
        }
        
        
        
        result_df <- data.frame(Comparison = pairs, adonis_pval = p_vals, adonis_R2 = r2s, p_val_adj = p.adjust(p_vals, p.adj.method), F.Model = F.Models)
        # --
        
        # - add vegan::adonis including all group levels -
        fit <- vegan::adonis(dist_obj ~ group_fac, permutations = nperm)
        df <- data.frame(Comparison = "Overall",
                         adonis_pval = fit$aov.tab[1, "Pr(>F)"], 
                         adonis_R2 = fit$aov.tab[1, "R2"], 
                         p_val_adj = fit$aov.tab[1, "Pr(>F)"], 
                         F.Model = fit$aov.tab[1, "F.Model"])
        
        result_df <- rbind(df, result_df)
        # -- 
        
        # - add significance levels -
        symnum.args$x <- result_df$adonis_pval
        result_df$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        result_df <- mutate(result_df, adonis_R2_PC = round(100*adonis_R2, 2))
        # --
        
}
# --







# --
#######################################
### FUNCTION: calc_ordination_from_distances
#######################################
# fully based on phyloseq::ordinate and phyloseq::plot_ordination


calc_ordination_from_distances <- function(physeq, group_var, dist_list, color_levels, ordination_type = "PCoA", shape = NULL, coord_cor = FALSE, phylum_colors = NULL, paired_var = NULL){
        
        
        
        # - prune physeq based on names(color_levels) and test that the inputs fit to each other -
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the phyloseq object.")
        }
        
        
        if(!is.null(paired_var) && !paired_var %in% colnames(sample_data(physeq))) {
                stop("The given paired_var is not a variable in the sample data of the phyloseq object.")
        }
        
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all names in names(color_levels)are actually levels in the group_var column.")
        }
        
        keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
        physeq <- prune_samples(keepSamples, physeq)
        
        if (!("dist" %in% class(dist_list[[1]]))){
                stop("first object of dist_list is not of class dist")
        }
        
        
        if (length(dist_list[[1]]) != nsamples(physeq)*(nsamples(physeq)-1)/2) {
                stop("The number of samples in the pruned physeq does not fit to the number of distances in dist_obj. Check and come back.")
        }
        
        # --
        
        ordination_list <- vector("list", length(dist_list))
        DFList <- vector("list", length(dist_list))
        DF_taxa_List <- vector("list", length(dist_list))
        # TrList <- vector("list", length(dist_list))
        TrList_own <- vector("list", length(dist_list))
        TrList_taxa <- vector("list", length(dist_list))
        
        axes <- 1:2 # currently only allowed to plot first and second
        
        for (i in seq_along(dist_list)) {
                
                ordination <- phyloseq::ordinate(physeq, method = ordination_type, distance = dist_list[[i]])
                ordination_list[[i]] <- ordination
                DF <- phyloseq::plot_ordination(physeq, ordination_list[[i]], color = group_var, justDF = TRUE)
                DFList[[i]] <- DF # just the first two axes cbind to sample_data in physeq
                
                x = colnames(DF)[1]
                y = colnames(DF)[2]
                Tr <- ggplot(DF, aes_string(x = x, y = y, col = group_var, shape = shape)) 
                Tr <- Tr + geom_point() + 
                        scale_color_manual("", values = color_levels) +
                        theme_bw() +
                        ggtitle(names(dist_list)[i])
                
                if (!is.null(paired_var)) {
                        Tr <- Tr + geom_line(aes_string(group = paired_var), col = cbPalette[1])
                }
                
                # for labelling axes
                if (ordination_type == "PCoA" || ordination_type == "NMDS") {
                        if (ordination_type == "PCoA") {
                                eigvec <- phyloseq:::extract_eigenvalue.pcoa(ordination)
                        } else {
                                eigvec <- phyloseq:::extract_eigenvalue.default(ordination)
                        }
                        
                        if (length(eigvec[axes]) > 0){
                                fracvar = eigvec[axes]/sum(eigvec)
                                percvar = round(100 * fracvar, 1)
                                strivar = as(c(Tr$label$x, Tr$label$y), "character")
                                strivar = paste0(strivar, " (", percvar, " %)")
                                Tr <- Tr + xlab(strivar[1]) + ylab(strivar[2]) 
                        }
                        
                        if (!is.null(eigvec) && coord_cor) {
                                Tr <- Tr + coord_fixed(sqrt(eigvec[2] / eigvec[1]))
                        }
                        
                }
                
                TrList_own[[i]] <- Tr
                rm(Tr)
                
                
                
                # TrList[[i]] <- phyloseq::plot_ordination(physeq, ordination_list[[i]], color = group_var) + ggtitle(names(dist_list)[i])
                
                DF_taxa <- phyloseq::plot_ordination(physeq, ordination_list[[i]], type = "taxa", color = "Phylum", justDF = TRUE)
                
                # - define or use phylum colors -
                if (is.null(phylum_colors)){
                        if (length(unique(DF_taxa$Phylum)) <= 15) {
                                fill_colors <- make_color_vector(DF_taxa$Phylum, rev(QuantColors15[1:length(unique(DF_taxa$Phylum))]))
                        } else {
                                fill_colors <- make_color_vector(DF_taxa$Phylum, viridis(length(unique(DF_taxa$Phylum))))
                        }
                        
                        DF_taxa$Phylum[is.na(DF_taxa$Phylum)] <- "NA" #pools all NA
                        DF_taxa$Phylum <- factor(DF_taxa$Phylum, levels = names(fill_colors), ordered = TRUE)
                        
                } else {
                        fill_colors <- phylum_colors
                        DF_taxa$Phylum[is.na(DF_taxa$Phylum)] <- "NA" #pools all NA
                        DF_taxa$Phylum <- factor(DF_taxa$Phylum, levels = names(fill_colors), ordered = TRUE)
                }
                
                # --
                
                
                DF_taxa_List[[i]] <- DF_taxa
                x = colnames(DF_taxa)[1]
                y = colnames(DF_taxa)[2]
                Tr <- ggplot(DF_taxa, aes_string(x = x, y = y, col = "Phylum")) 
                Tr <- Tr + geom_point() +
                        scale_color_manual("", values = phylum_colors) +
                        theme_bw() +
                        ggtitle(names(dist_list)[i])
                
                # for labelling axes
                if (ordination_type == "PCoA" || ordination_type == "NMDS") {
                        if (ordination_type == "PCoA") {
                                eigvec <- phyloseq:::extract_eigenvalue.pcoa(ordination)
                        } else {
                                eigvec <- phyloseq:::extract_eigenvalue.default(ordination)
                        }
                        
                        if (length(eigvec[axes]) > 0){
                                fracvar = eigvec[axes]/sum(eigvec)
                                percvar = round(100 * fracvar, 1)
                                strivar = as(c(Tr$label$x, Tr$label$y), "character")
                                strivar = paste0(strivar, " (", percvar, " %)")
                                Tr <- Tr + xlab(strivar[1]) + ylab(strivar[2]) 
                        }
                        
                        if (!is.null(eigvec) && coord_cor) {
                                Tr <- Tr + coord_fixed(sqrt(eigvec[2] / eigvec[1]))
                        }
                        
                }
                
                TrList_taxa[[i]] <- Tr
                rm(Tr)
                
        }
        
        names(ordination_list) <- names(TrList_taxa) <- names(DFList) <- names(DF_taxa_List) <- names(TrList_own) <- names(dist_list)
        out <- list(ordination_list = ordination_list, DFList = DFList, DF_taxa_List = DF_taxa_List, ordination_Tr_samples = TrList_own, ordination_Tr_taxa = TrList_taxa)
}
# --





# --
#######################################
### FUNCTION: test_diffs_in_prevalence_single
#######################################
# Function performs fisher exact test on prevalence (absence presence) between the levels
# in a grouping factor
# INPUT:
# physeq: phyloseq
# group_var: name of the column in sample_data(physeq) that defines the groups
# p.adj.method, used in p.adjust
# minCount: present are taxa in species with more counts than minCount
# OUTPUT:
# list of data.frames, one data frame for each combi of levels in your grouping factor
# The data frames are ordered by p_value, and the tax_table has been cbound:)

test_diffs_in_prevalence_single <- function(physeq, group_var, compare = NULL, p.adj.method = "fdr", minCount = 0L, symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        
        CT <- as(otu_table(physeq), "matrix")
        CT <- CT > minCount
        
        
        prev_list <- lapply(group_var_levels, function(level){
                data.frame(Present = colSums(CT[group_fac == level, ]),
                           Absent = colSums(!CT[group_fac == level, ]))
        })
        
        pr_ab_gr1 <- prev_list[[1]] #pr_ab = presence absence
        pr_ab_gr2 <- prev_list[[2]]
        
        rowwise_compare_matrix <- cbind(pr_ab_gr1, pr_ab_gr2)
        
        FisherTests <- lapply(1:nrow(rowwise_compare_matrix), function(e){
                mat_fisher <- matrix(c(rowwise_compare_matrix[e, 1],
                                       rowwise_compare_matrix[e, 3],
                                       rowwise_compare_matrix[e, 2],
                                       rowwise_compare_matrix[e, 4]), ncol = 2)
                fisher.test(mat_fisher, conf.int = TRUE, conf.level = 0.95)
                # fisher.test(mat_fisher, conf.int = FALSE)
        })
        
        # - take out oddsRatios -
        oddRatios <- sapply(FisherTests, function(TestResult){TestResult$estimate})
        oddRatios_lb <- sapply(FisherTests, function(TestResult){TestResult$conf.int[1]})
        oddRatios_ub <- sapply(FisherTests, function(TestResult){TestResult$conf.int[2]})
        
        direction <- rep(group_var_levels[1], length(oddRatios))
        direction[oddRatios < 1] <- group_var_levels[2]
        
        oddRatios_lb[oddRatios < 1] <- 1/oddRatios_lb[oddRatios < 1]
        oddRatios_ub[oddRatios < 1] <- 1/oddRatios_ub[oddRatios < 1]
        oddRatios_lb_final <- pmin(oddRatios_lb, oddRatios_ub)
        oddRatios_ub_final <- pmax(oddRatios_lb, oddRatios_ub)
        oddRatios[oddRatios < 1] <- 1/oddRatios[oddRatios < 1]
        # --
        
        
        # - take out p_vals and assign significance levels -
        p_vals <- sapply(FisherTests, function(TestResult){TestResult$p.value})
        p_vals_adj <- p.adjust(p_vals, p.adj.method)
        
        symnum.args$x <- p_vals
        significance <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- p_vals_adj
        significance_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
        # --
        
        # - take out prevalence percentages -
        prev_PC_gr1 <- round(100*(pr_ab_gr1[, "Present"]/sum(group_fac == group_var_levels[1])), 1)
        prev_PC_gr2 <- round(100*(pr_ab_gr2[, "Present"]/sum(group_fac == group_var_levels[2])), 1)
        # --
        
        df <- data.frame(p_val = p_vals, p_val_adj = p_vals_adj, signi = significance, signi_adj = significance_adj, oddsRatio = round(oddRatios, 2), oddsRatio_lb = round(oddRatios_lb_final, 2),
                         oddsRatio_ub = round(oddRatios_ub_final, 2), direction = direction, comparison = paste(group_var_levels, collapse = " vs "),  prev_PC_gr1 = prev_PC_gr1,  prev_PC_gr2 = prev_PC_gr2)
        
        df <- cbind(as.data.frame(df), as.data.frame(unclass(tax_table(physeq))))
        df$Taxon <- colnames(CT)
        df <- arrange(df, p_val) %>% select(Taxon, 1:(ncol(df)-1))
        df
        
}
# --






# --
#######################################
### FUNCTION: format_hit_table
#######################################

format_hit_table <- function (result_df, p.adjust.threshold = 0.1, p.adjust.method = NULL) {
        
        if (!all(c("p_val", "p_val_adj", "Taxon", "direction", "signi", "signi_adj") %in% colnames(result_df))) {
                stop("result_df should be a data.frame generated by one of the differential abundance tests and contain all corresponding columns.")
        }
        
        if (!is.null(p.adjust.method)) {
                result_df$p_val_adj = p.adjust(result_df$p_val, method = p.adjust.method)
        }
        
        
        result_df <- arrange(result_df, p_val_adj, p_val)
        
        no_hits <- sum(result_df$p_val_adj <= p.adjust.threshold, na.rm = TRUE)
        
        keepTaxa <- no_hits
        
        if (keepTaxa < 10 && nrow(result_df) >= 10) {
                keepTaxa <- 10
        } else if (keepTaxa < 10 && nrow(result_df < 10)){
                keepTaxa <- nrow(result_df)
        }
        
        
        df <- result_df[1:keepTaxa,]
        
        taxa_annotation <- get_taxon_names(df)
        taxa_annotation <- strsplit(taxa_annotation, split = "/")
        taxa_annotation <- sapply(taxa_annotation, `[`, 1)
        df$Annotation <- taxa_annotation
        
        df <- select(df, Taxon, Annotation, p_val, p_val_adj, signi, signi_adj, direction, colnames(df)[!(colnames(df) %in% c("Taxon", "Annotation", "p_val", "p_val_adj", "signi", "signi_adj", "direction"))])
        
        rownames(df) <- df$Taxon
        
        list(hit_table = df, no_hits = no_hits)
}
# --


# --
#######################################
### FUNCTION: get_taxon_names
#######################################
# NB: adjusted it here for CRC course to only show the actual level

get_taxon_names <- function(df) {
        
        df1 <- df[, colnames(df) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", 
                                        "Genus", "Species")]
        
        if (ncol(df1) == 0) {stop("the provided data frame did not contain the expected taxonomy columns such as Phylum, Class etc.")}
        
        df1[] <- lapply(df1, as.character)
        
        Last_NotNA_Position <- apply(df1, 1, function(x){length(which(!is.na(x)))})
        Last_NotNA_Position[Last_NotNA_Position == 0] <- 1
        
        Names <- vector(mode = "character", length = nrow(df1))
        
        for (i in 1:nrow(df1)) {
                #if (Last_NotNA_Position[i] == 7){
                 #       Names[i] <- paste(df1[i, "Genus"], df1[i, "Species"], sep = " ")
                #} else {
                        Names[i] <- df1[i, Last_NotNA_Position[i]] 
                #}
        }
        
        Names[is.na(Names)] <- "NA"
        Names
        
}
# --





# --
#######################################
### test_differential_abundance_DESeq2single
#######################################
## Inputs
# physeq: phyloseq object
# group_var: name of column that defines group fac in sample_data
# SFs: often you might want to give the SizeFactors already because you wanted to calculate them on non-filtered data,
# when SFs are not NULL, type is ignored
# type: type in estimateSizeFactors, ignored when Size factors given
## OUTPUT:
# list, first: result df of DESEQ2 analysis, second: the adjusted phyloseq object after size factor correction


# ATTENTION: you could add here block as Mani had in test_differential_abundance_DESeq2

test_differential_abundance_DESeq2single <- function(physeq, group_var, compare = NULL, cooksCutoff = TRUE, SFs = NULL, type = "ratio", p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))){
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        
        # NB: DESeq2 can not deal with ordered factors, sees them somehow as just one level, therefore
        sample_data(physeq)[[group_var]] <- factor(group_fac, levels = c(group_var_levels, setdiff(levels(group_fac), group_var_levels)), ordered = FALSE)
        
        DES = phyloseq::phyloseq_to_deseq2(physeq, formula(paste("~", group_var)))
        
        
        if (is.null(SFs)){
                if(type == "ratio"){
                        GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = FALSE)
                }
                
                dds <- estimateSizeFactors(DES, type = type, geoMeans = GM)
                # NB: geoMeans is ignored when type = "iterate"
                # NB2: "iterate" takes much longer than "ratio", and when using GM via gm_own with "ratio" the size factors 
                # correlate with above 99% with size factors from "iterate"
                # NB3: "poscounts" is the same as "ratio" with GM calculated with zeros.count = TRUE!!
                # SFs2 <- sizeFactors(dds)
                
        } else {
                dds <- DES
                sizeFactors(dds) = SFs
                # identical(sizeFactors(dds), SFs) 
        }
        
        
        dds <-  DESeq(dds, fitType = "parametric", test = "Wald", 
                      quiet = TRUE, minReplicatesForReplace = Inf) 
        
        # to get the size factor adjusted physeq object
        physeq_out <- physeq
        otu_table(physeq_out) <- otu_table(t(counts(dds, normalized = TRUE)), taxa_are_rows = FALSE)
        
        
        
        # - analyse the res -
        res <- as.data.frame(results(dds, contrast = c(group_var, group_var_levels), cooksCutoff = cooksCutoff))
        
        res$p_val_adj <- p.adjust(res$pvalue, method = p.adjust.method) # NB: in case of "fdr" same as default DESeq2
        
        CT <- counts(dds, normalized = TRUE)
        i = 1
        j = 2
        n1 <- sum(group_fac == group_var_levels[i])
        n2 <- sum(group_fac == group_var_levels[j])
        # res$Median_grp1 <- apply(CT[, group_fac == group_var_levels[i]], 1, median)
        # res$Median_grp2 <- apply(CT[, group_fac == group_var_levels[j]], 1, median)
        res$Mean_grp1 <- apply(CT[, group_fac == group_var_levels[i]], 1, mean)
        res$Mean_grp2 <- apply(CT[, group_fac == group_var_levels[j]], 1, mean)
        # res$baseMeanSelf <- apply(CT, 1, mean) # exactly the same as baseMean!
        # res$Zeros_grp1 <- apply(CT[, group_fac == group_var_levels[i]], 1, function(cnts){sum(cnts == 0)})
        # res$Zeros_grp2 <- apply(CT[, group_fac == group_var_levels[j]], 1, function(cnts){sum(cnts == 0)})
        res$Present_grp1 <- apply(CT[, group_fac == group_var_levels[i]], 1, function(cnts){sum(cnts != 0)})
        res$Present_grp2 <- apply(CT[, group_fac == group_var_levels[j]], 1, function(cnts){sum(cnts != 0)})
        res$prev_PC_grp1 <- round(100*(res$Present_grp1/n1),1)
        res$prev_PC_grp2 <- round(100*(res$Present_grp2/n2), 1)
        #res$Sparsity_grp1 <- 100*(res$Zeros_grp1/n1)
        #res$Sparsity_grp2 <- 100*(res$Zeros_grp2/n2)
        
        # - add sginificance and direction -
        symnum.args$x <- res$pvalue
        res$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- res$p_val_adj
        res$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character() # note ? in case of p_val = NA
        
        res$direction <- rep(group_var_levels[2], nrow(res))
        res$direction[res$log2FoldChange > 0] <- group_var_levels[1]
        
        res$comparison <- paste(group_var_levels, collapse = " vs ")
        # --
        
        res$Taxon <- rownames(res)
        
        res <- dplyr::select(res, Taxon, teststat = stat, p_val = pvalue, p_val_adj,
                             signi, signi_adj, direction, comparison,
                             baseMean, log2FoldChange, Mean_grp1,
                             Mean_grp2, prev_PC_grp1, prev_PC_grp2)
        
        # NB: I dropped here padj from DESeq since same as p_val_adj in case of p.adjust.method = "fdr"
        res <- cbind(res, as.data.frame(unclass(tax_table(physeq))))
        res <- dplyr::arrange(res, desc(abs(teststat)))
        list(result_df = res, physeq_out = physeq_out) 
        
}
# --



# --
#######################################
### test_differential_abundance_Wilcoxonsingle ##
#################

test_differential_abundance_Wilcoxonsingle <- function(physeq, group_var, compare = NULL, excludeZeros = FALSE, p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))){
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        CT <- as(otu_table(physeq), "matrix")
        
        i <- group_var_levels[1]
        j <- group_var_levels[2]
        
        res_mat <- apply(CT, 2, function(taxon_counts){
                x <- taxon_counts[group_fac == i]
                Zeros_grp1 <- sum(x == 0)
                # # Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                Present_grp1 <- length(x)-Zeros_grp1
                prev_PC_grp1 <- 100*(Present_grp1/length(x))
                if(excludeZeros){
                        x <- x[x != 0]
                }
                Median_grp1 <- median(x, na.rm = T) # NA in case all 0 and excludeZeros was TRUE
                Mean_grp1 <- mean(x, na.rm = T) # NaN in case all 0 and excludeZeros was TRUE
                if (is.na(Mean_grp1)){ Mean_grp1 = NA }
                
                y <- taxon_counts[group_fac == j]
                Zeros_grp2 <- sum(y == 0)
                Present_grp2 <- length(y)-Zeros_grp2
                # # Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                prev_PC_grp2 <- 100*(Present_grp2/length(y))
                if(excludeZeros){
                        y <- y[y != 0]
                }
                Median_grp2 <- median(y, na.rm = T)
                Mean_grp2 <- mean(y, na.rm = T)
                if (is.na(Mean_grp2)){ Mean_grp2 = NA }
                
                if (length(x) != 0 && length(y) != 0){
                        wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                        pValue <- wilcTest$p.value
                        W <- wilcTest$statistic
                        # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                        Ranks <- rank(c(x, y))
                        n1 <- length(x)
                        n2 <- length(y)
                        # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                        # how about the other W?
                        # Wy <- sum(Ranks[(n1+1):(n1+ n2)]) - (n2*(n2+1)/2)
                        standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                        
                        # # if you want to check that multtest::mt.minP would give the same statistic
                        # mati <- matrix(c(x,y), nrow = 1)
                        # grFac <- c(rep(fac_levels[i], n1), rep(fac_levels[j], n2))
                        # grFac <- factor(grFac, levels = c(fac_levels[i], fac_levels[j]))
                        # standStat2 <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                        # # identical(standStat, standStat2) # TRUE
                        # uncomment all with standStat2 to test all the way
                        
                } else {
                        pValue = NA
                        W <- NA
                        standStat = NA
                        n1 <- length(x)
                        n2 <- length(y)
                        # standStat2 = NA
                }
                
                
                c(standStat, pValue, Median_grp1, Median_grp2, 
                  Mean_grp1, Mean_grp2, prev_PC_grp1, prev_PC_grp2, n1, n2, W) 
        })
        
        res_mat <- t(res_mat)
        
        DF <- data.frame(Taxon = rownames(res_mat), res_mat)
        colnames(DF) <- c("Taxon", "teststat", "p_val", "Median_grp1", "Median_grp2", "Mean_grp1", "Mean_grp2", "prev_PC_grp1", "prev_PC_grp2", "n1", "n2", "W")
        DF$p_val_adj <- p.adjust(DF$p_val, method = p.adjust.method)
        
        
        symnum.args$x <- DF$p_val
        DF$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- DF$p_val_adj
        DF$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
        
        DF$direction <- i
        DF$direction[DF$teststat > 0] <- j
        DF$comparison <- paste(group_var_levels, collapse = " vs ")
        
        
        DF <- dplyr::select(DF,  Taxon:p_val, p_val_adj:comparison, Median_grp1:W)
        
        DF <- cbind(DF, as.data.frame(unclass(tax_table(physeq))))
        DF <- dplyr::arrange(DF, desc(abs(teststat)))
        DF
}
# --






# --
####################################
## plot_taxa_ratios_AllLevels 
###################################
# see plot_taxa_ratios_levelPairs: Here you directly calculate the count by count ratio matrix only for the tax_nom, you still facet by taxa_den
# (denominator) but you keep all levels in all plots. Makes extensive use of ggpubr, NB: ggpubr is so smart to adjust p-values when you use
# scale_y_log10

## Output: 
# - list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)


plot_taxa_ratios_AllLevels <- function(physeq, group_var, color_levels, tax_names = NULL,
                                       taxa_nom = "Firmicutes", taxa_den = NULL, test = "t.test", p_adjust_method = "fdr",
                                       tax_order = NULL,
                                       symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), hide.ns = FALSE) {
        
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all names in names(color_levels)are actually levels in the group_var column.")
        }
        
        if (!(test %in% c("t.test", "wilcox.test"))) {
                stop("test should be t.test or wilcox.test")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        compare <- names(color_levels)
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        # if (length(group_var_levels) != 2) {
        #         stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
        #                     paste(group_var_levels, collapse = ", ")))
        # }
        
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:ntaxa(physeq), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        taxa_names(physeq) <- tax_names
        # --
        
        # - calculate the matrix taxa_nom/(all other taxa) -         
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        # ONLY keep samples defined by group_var_levels (= names(color_levels))
        CT <- CT[, group_fac %in% group_var_levels]
        
        
        i <- which(rownames(CT) == taxa_nom)
        if (length(i) != 1) {stop("taxa_nom not found in tax_names or tax_names not unique!")}
        
        
        TbTmatrix <- apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are possibly Inf, and NaN values in the matrix, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        # --
        
        TbT_DF <- as.data.frame(TbTmatrix)
        TbT_DF$Taxon <- rownames(TbT_DF)
        # - use taxa_den (denominator) to restrict the taxa to which taxa_nom is compared to -
        if (is.null(taxa_den)) {taxa_den <- tax_names}
        TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
        # --
        
        # - change to long DF -
        TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
        # --
        
        # - add the group_var level information  -
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        TbT_DF_l$Group <- as.character(LookUpDF$Group[match(TbT_DF_l$Sample, LookUpDF$Sample)])
        TbT_DF_l$Group <- factor(TbT_DF_l$Group, levels = group_var_levels, ordered = T)
        # --
        
        # - change all ratios where either nominator taxon or denominator taxon had count = 0 to NA -
        # remember: 0/0 = NaN (not finite), 0/x = 0, x/0 = Inf
        TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
        # --
        
        # - find taxa that would throw an error in statistical test and remove those taxa from DF -
        # first find the taxa that would throw an error in t.test or wilcox.test
        var_plus_length_check <- group_by(TbT_DF_l, Taxon, Group) %>% summarise(Variance = var(Ratio, na.rm = T), NotNA = sum(!is.na(Ratio)))
        if (test == "t.test"){
                var_plus_length_check <- dplyr::filter(var_plus_length_check, !(Variance > 0) | NotNA < 2) # variance > 0 also to remove test where host_taxon == taxon
        } else if (test == "wilcox.test") {
                var_plus_length_check <- dplyr::filter(var_plus_length_check, !(Variance > 0) | NotNA < 1)
        }
        
        if (nrow(var_plus_length_check) != 0) { 
                TbT_DF_l <- filter(TbT_DF_l, !(Taxon %in% unique(var_plus_length_check$Taxon)))
        }
        # --
        
        # - use ggpubr::compare_means to calculate all pValues of the Ratios for the different taxa_den between current group_levels -
        pVals <- ggpubr::compare_means(formula = Ratio ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        
        pVals <- dplyr::arrange(pVals, p)
        
        # NB: the pVals for t.test change when you log the ratios (scale_y_log10())
        TbT_DF_l$RatioLog10 <- log10(TbT_DF_l$Ratio)
        pValsLog <- ggpubr::compare_means(formula = RatioLog10 ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        pValsLog <- dplyr::arrange(pValsLog, p)
        # --
        
        # NB: I plot now for log and non log independently, even though ggpubr is so smart to change the p-values when you just use
        # Tr + scale_y_log10(). I plot independently because log might change the order in case of t.test!
        TbT_DF_l_log <- TbT_DF_l
        
        # - order taxa based on pVals result or based on tax_order -
        if (is.null(tax_order)){
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = unique(pVals$Taxon), ordered = TRUE)
                TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = unique(pValsLog$Taxon), ordered = TRUE)
                
        } else {
                if(!all(unique(pVals$Taxon) %in% tax_order)){
                        stop("given tax_order does not fit to tax_names")
                }
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = tax_order, ordered = TRUE)
                TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = tax_order, ordered = TRUE)
        }
        # --
        
        # - since you might have more than two levels in each plot you need to set the comparisons argument in stat_compare_means -
        comparisonList <- get_unique_facLevel_combinations(group_var_levels)
        # --
        
        # - plot: NB: also for non removed taxa some samples might have NA ratios that will be removed -
        
        # Tr <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
        # Tr <- Tr +
        #         geom_violin() +
        #         geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
        #         # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
        #         scale_color_manual(values = color_levels) +
        #         facet_wrap(~ Taxon, scales = "free_y") +
        #         xlab("") +
        #         ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
        #         theme_bw() +
        #         theme(legend.position = "none")
        # 
        # Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        # # Tr <- Tr + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5, hide.ns = hide.ns)
        # 
        # Tr1 <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
        # Tr1 <- Tr1 +
        #         geom_boxplot(outlier.color = NA) +
        #         geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
        #         # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
        #         scale_color_manual(values = color_levels) +
        #         facet_wrap(~ Taxon, scales = "free_y") +
        #         xlab("") +
        #         ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
        #         theme_bw() +
        #         theme(legend.position = "none")
        # 
        # Tr1 <- Tr1 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        # #Tr1 <- Tr1 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns) # p.format
        
        
        # Tr2 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
        # Tr2 <- Tr2 +
        #         geom_violin() +
        #         geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
        #         # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
        #         scale_color_manual(values = color_levels) +
        #         facet_wrap(~ Taxon, scales = "free_y") +
        #         xlab("") +
        #         ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
        #         theme_bw() +
        #         theme(legend.position = "none")
        # 
        # 
        # 
        # # Tr2 <- Tr2 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns)
        # Tr2 <- Tr2 + scale_y_log10()
        # Tr2 <- Tr2 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        
        
        Tr3 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
        Tr3 <- Tr3 +
                geom_boxplot(outlier.color = NA) +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("Abundance ratio: ", taxa_nom, "/stated phylum", sep = "")) +
                theme_bw() +
                theme(legend.position = "none")
        
        # Tr3 <- Tr3 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns)
        Tr3 <- Tr3 + scale_y_log10()
        Tr3 <- Tr3 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        
        
        # list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)
        
        list(pVals = pVals, pValsLog = pValsLog, Tr3 = Tr3)
        
}
#--




# --
####################################
## calculate_raw_TbTmatrixes:
###################################

calculate_raw_TbTmatrixes = function(physeq){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        
        TbTmatrixes <- lapply(1:nrow(CT), function(i){apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are Inf, and NaN values in the matrixes, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        
        names(TbTmatrixes) <- rownames(TbTmatrixes[[1]])
        
        
        TbTmatrixes
        
}
# --








# --
####################################
## create_raw_TbT_TilePlot: 
###################################
# NB: In this version TbTmatrixes should have been calculated on all samples in physeq. 
# names(color_levels) should only contain two levels in group_var!

create_raw_TbT_TilePlot <- function(TbTmatrixes, physeq, group_var, color_levels, tax_names = NULL, tax_order = NULL, 
                                    test = "wilcoxon", signi_level = 0.05, p_adjust_method = "none") {
        
        if(!identical(length(TbTmatrixes), ntaxa(physeq))){stop("TbTmatrixes don't fit to physeq")}
        
        if(ncol(TbTmatrixes[[1]]) != nsamples(physeq)){stop("TbTmatrixes don't fit to physeq. Not the same number of samples.")}
        
        
        if(test != "wilcoxon" & test != "t.test"){stop("test unknown, must be wilcoxon or t.test")}
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        compare <- names(color_levels)
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (names(color_levels)) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:ntaxa(physeq), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        # --
        
        
        
        names(TbTmatrixes) <- tax_names
        
        TbTmatrixes <- lapply(TbTmatrixes, function(mat){
                rownames(mat) <- tax_names
                mat
        })
        
        
        i <- group_var_levels[1]
        j <- group_var_levels[2]
        
        # ntaxa * ntaxa wilcoxon tests take time if you have a lot of taxa!
        pValMatrix <- sapply(TbTmatrixes, function(mat){
                apply(mat, 1, function(taxon_ratios){
                        x <- taxon_ratios[group_fac == i]
                        x <- x[is.finite(x) & x != 0] # removes all ratios in which one of the two taxa was not present!
                        y <- taxon_ratios[group_fac == j]
                        y <- y[is.finite(y) & y != 0] # removes 0/0 = NaN, 0/x = 0, x/0 = Inf
                        if (test == "wilcoxon"){
                                if (length(x) > 0 && length(y) > 0){
                                        pValue <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)$p.value
                                        # NB: wilcox.test ignores in default setting (maybe see na.action) NA, NaN, Inf, -Inf
                                        # For plot: change sign of pValue to negative if taxon is more abundant in group 1. 
                                        Ranks <- rank(c(x[!is.na(x)], y[!is.na(y)]))
                                        n1 <- length(x[!is.na(x)])
                                        n2 <- length(y[!is.na(y)])
                                        Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2)
                                        Wy <- sum(Ranks[(n1+1):(n1+n2)])-(n2*(n2+1)/2)
                                        if(Wx > Wy){pValue <- -1*pValue}
                                        pValue
                                        
                                } else {
                                        pValue = 1
                                }
                                
                        } else if (test == "t.test") {
                                if (length(x) > 1 && length(y) > 1 && var(x) > 0 && var(y) > 0){
                                        pValue <- t.test(x = x, y = y, alternative = "two")$p.value
                                        if (mean(x, na.rm = T) > mean(y, na.rm = T)){pValue <- -1*pValue}
                                        pValue
                                } else {
                                        pValue <- 1
                                }
                                
                        }
                        
                })
        })
        
        # make sure diagonal is all NA (can be exceptions especially for t.test)
        diag(pValMatrix) <- NA
        
        # - adjust p-values if asked for -        
        signs <- pValMatrix < 0
        signs[is.na(signs)] <- FALSE
        
        pValMatrix <- abs(pValMatrix)
        for (e in 1:nrow(pValMatrix)){
                pValMatrix[e, ] <- p.adjust(pValMatrix[e, ], method = p_adjust_method)
        } # equal to t(apply(pValMatrix, 1, p.adjust, method = p_adjust))
        
        pValMatrix[signs] <- pValMatrix[signs]*(-1)
        # --
        
        # -- add a tile plot of the pValMatrix --
        
        DF <- as.data.frame(pValMatrix)
        DF[is.na(DF)] <- 2 # just to avoid missing values in plot and have a clear non-pValue value to mark self comparisons as black
        DF$HostTaxon <- rownames(pValMatrix)
        DF <- tidyr::gather(DF, key = Taxon , value = pValue, - HostTaxon)
        if (is.null(tax_order)) {
                DF$Taxon <- factor(DF$Taxon, levels = rownames(pValMatrix), ordered = TRUE)
                DF$HostTaxon <- factor(DF$HostTaxon, levels = rev(rownames(pValMatrix)), ordered = TRUE)
        } else {
                if(!all(rownames(pValMatrix) %in% tax_order)){
                        stop("given tax_order does not fit to tax_names")
                }
                DF$Taxon <- factor(DF$Taxon, levels = tax_order, ordered = TRUE)
                DF$HostTaxon <- factor(DF$HostTaxon, levels = rev(tax_order), ordered = TRUE)
                
        }
        
        fill_colors <- c(color_levels, ns = "gray98", "test not possible" = "black")
        
        DF$Fill <- "ns"
        DF$Fill[DF$pValue < signi_level & DF$pValue > 0] <- i
        DF$Fill[DF$pValue > -1*signi_level & DF$pValue < 0] <- j
        DF$Fill[DF$pValue == 2] <- "test not possible"
        DF$Fill <- factor(DF$Fill, levels = names(fill_colors), ordered = T)
        TileTr <- ggplot(DF, aes(x = Taxon, y = HostTaxon, fill = Fill))
        TileTr <- TileTr + 
                geom_raster() + 
                # ggtitle(paste(i, " vs ", j, sep = "")) +
                scale_fill_manual("", values = fill_colors) +
                scale_x_discrete(position = "top") +
                labs(x=NULL, y=NULL) +
                theme_bw() +
                #theme_tufte(base_family="Helvetica") +
                theme(panel.border = element_blank(),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
                      axis.ticks=element_blank())
        
        TileTr
        
        
}
# --







