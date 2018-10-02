# How I changed the microdiab phyloseq object

datapath <- "/Users/jvb740/Coursera_MOOC/20161202_LearningShiny_FantasySports/shinyy/Apps/Teaching_Apps/MicrobiomeX2_App"

phyloseq_name <- "physeq_microdiab_ngt_Men.rds"


ps <- readRDS(file = file.path(datapath, phyloseq_name))

# - added to remove taxa that are not present in a single sample directly -
keepTaxa <- taxa_names(ps)[taxa_sums(ps) > 0]
ps <- phyloseq::prune_taxa(keepTaxa, ps)
# --


# - change sample_data() -
SS <- as(sample_data(ps), "data.frame")
SS$Country <- as.character(SS$Country)
SS$Country[SS$Country == "IN"] <- "SL"
SS$ID <- paste0("Sample_", 1:nrow(SS))
set.seed(12334)
SS$Gender <- sample(x = c("M", "F"), size = nrow(SS), replace = T)
SS$Status <- "Healthy"
SS$Age <- sample(x = 19:72, size = nrow(SS), replace = T)
SS <- dplyr::select(SS, ID:Status, Age)
rownames(SS) <- SS$ID
# --

# - change tax-table -
TT <- as.data.frame(as(tax_table(ps), "matrix"))
rownames(TT) <- paste0("T_", sprintf('%0.4d', 1:nrow(TT)))

TT[] <- lapply(TT, as.character)

# checked prevalence to find 2 that had ca 40% prevalence
TT[901,6] <- "Fakides"
TT[901,7] <- "legitimus"
TT[60,6] <- "Fakides"
TT[60,7] <- "fatalis"
# --

# - change otu_table -
OTU <- as(otu_table(ps), "matrix")
colnames(OTU) <- rownames(TT)
rownames(OTU) <- rownames(SS)
EffectSize <- 150
OTU[SS$Country == "DK",901] <- OTU[SS$Country == "DK",901] * EffectSize
OTU[SS$Country == "SL",60] <- OTU[SS$Country == "SL",60] * EffectSize
# --


# - generate phyloseq object -
OTU <- otu_table(OTU, taxa_are_rows = FALSE)
TT <- as(TT, "matrix")
TT <- tax_table(TT)
SS <- sample_data(SS)

ps_new <- phyloseq(OTU, SS, TT)
# --

saveRDS(object = ps_new, "Project_PSobject_CountryComparison.rds")





