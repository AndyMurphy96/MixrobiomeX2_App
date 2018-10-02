# source("https://bioconductor.org/biocLite.R")
# first time these bioconductor packages have to be installed with biocLite("phyloseq")
library(phyloseq)
library(vegan)
library(DESeq2)

# packages from cran have to be installed with e.g. install.packages("ggplot2")
library(ggplot2)
library(shiny)
library(lubridate)
library(tidyr)
library(dplyr)
library(viridis)
# library(rdrop2)
#library(xtable)
library(pheatmap)
# library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)

functionpath <- "./Functions"
#source(file.path(functionpath, "_nN_000_helper_functions.R"))
#source(file.path(functionpath, "TeachingFunctions.R"))
source(file.path(functionpath, "_n_000_helper_functions.R"))
source(file.path(functionpath, "_n_010_explore_ps_functions.R"))
source(file.path(functionpath, "_n_020_alpha_diversity_functions.R"))
source(file.path(functionpath, "_n_030_preprocess_filtering_functions.R"))
source(file.path(functionpath, "_n_040_beta_diversity_functions.R"))
source(file.path(functionpath, "_n_050_diff_abundance_functions.R"))
source(file.path(functionpath, "_n_060_phylum_analysis_functions.R"))
source(file.path(functionpath, "shiny_app_functions.R"))




ui <- fluidPage(
        tags$head(
                # set color of the info texts
                tags$style(HTML("
                                .shiny-text-output {
                                color: blue;
                                font-size: 15px;
                                }
                                ")),
                # set colors of checkboxes (not in currently)
                tags$style(HTML("
                                .checkbox {
                                color: black;
                                font-size: 22px;
                                line-height: 35px;
                                }
                                ")),
                # next set's text of the radio buttons
                tags$style(HTML("
                                .radio {
                                color: black;
                                font-size: 15px;
                                }
                                ")),
                # sets text of navbarPage title
                # rgb(17, 119, 85), red #cc3f3f
                tags$style(type = 'text/css',
                           '.navbar-default .navbar-brand {
                             color: #814681;
                                font-weight: 800;
                                font-size: 20px;
                           }',
                # sets background-color of navbarPage bar and text size of the tabs, not sure how to change the font color
                   '.navbar {
                        background-color: #23373B;
                        font-size: 15px;
                        font-weight: 700;
                }')


        ),
        
        
     
        
        #titlePanel("Microbiome Quiz"),
        navbarPage(title = "MicrobiomeX2",
                   
                   ############# Tab 1 ########################
                   
                   tabPanel(title = "Parameters",
                            sidebarLayout(
                                    sidebarPanel(
                                            
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Set general input parameters")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_1')
                                            ),
                                            
                                            wellPanel(
                                                    
                                                    tags$h4("Group variables and colors"),
                                                    textInput(inputId = "group_var", label = "group_var", value = "Country", width = "150px"),
                                                    splitLayout(
                                                            textInput(inputId = "grp1", label = "group 1", value = "DK", width = "120px"),
                                                            textInput(inputId = "grp2", label = "group 2", value = "SL", width = "120px")),
                                                    splitLayout(
                                                            textInput(inputId = "grp1_color", label = "group 1 color", value = "#009E73", width = "120px"),
                                                            textInput(inputId = "grp2_color", label = "group 2 color", value = "#D55E00", width = "120px")
                                                    )
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            tableOutput(outputId = "overviewView1")
                                            
                                            #tableOutput(outputId = "phylaViews"),
                                            
                                            #plotOutput(outputId = "samplesPS")
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 2 ########################
                   tabPanel(title = "Load and explore phyloseq object",
                            sidebarLayout(
                                    sidebarPanel(
                                            
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Get an overview of the data in the phyloseq object:"),
                                                    tags$ul(
                                                            tags$li("Step 1: Load phyloseq object"),
                                                            tags$li("Step 2: Get an idea of the phyla and their abundance distribution (in the raw data!). This step sets phylum_colors for later plots."),
                                                            tags$li("Step 3: Visualize sample sizes (library sizes/total counts).")
                                                    )
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_2')
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Step1: Load phyloseq object"),
                                                    
                                                    fileInput(inputId = "loadPhyseq", label = ""), # accept = "text/rds"
                                                    
                                                    tags$h5(),
                                                    wellPanel(
                                                            tags$h4("Look at the components of the phyloseq object"),
                                                            radioButtons(inputId = "component", label = "Select component of phyloseq object",
                                                                         choices = c("otu_table",
                                                                                     "sample_data",
                                                                                     "tax_table"),
                                                                         width = "100%"),
                                                            textInput(inputId = "sampleRange", label = "Sample range", value = "1:100"),
                                                            textInput(inputId = "taxaRange", label = "Taxa range", value = "1:100"),
                                                            actionButton(inputId = "checkPhyseq", label = "See component of phyloseq object")
                                                    )
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Select data to use for steps 2 and 3"),
                                                    
                                                    radioButtons(inputId = "tcaType2", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered2", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%")
                                                    
                                            ),
                        
                                            
                                            wellPanel(
                                                    tags$h4("Step 2: Phyla distribution"),
                                                 
                                                    
                                                    actionButton(inputId = "calcPhyla", label = "Calc. phyla distribution and phylum_colors"),
                                                    
                                                    tags$h5("NB: Step is necessary to calculate phylum_colors for later plots! Can be repeated after size factor adjustment.")
                                                    ),
                                            
                                            wellPanel(
                                                    tags$h4("Step 3: Overview of sample sizes"),
                                                    
                                                    actionButton(inputId = "barplotPS", label = "Generate abundance barplot of the samples"),
                                                    
                                                    tags$h5("NB: Please be patient, plot calculation takes 20 to 30 seconds.")
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            tableOutput(outputId = "overviewView2"),
                                            
                                            # tableOutput(outputId = "explorationViews"),
                                            
                                            tableOutput(outputId = "phylaViews"),
                                            
                                            plotOutput(outputId = "samplesPS")
                                            
                                            )
                                            # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 3 ########################
                   tabPanel(title = "Library size adjustment",
                            sidebarLayout(
                                    sidebarPanel(
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Adjust differences in library size of the samples, so samples can be compared to each other. Three variations of calculating size factors are implemented."),
                                                    tags$ul(
                                                            tags$li(tags$b("DESeq_gm_exclZero:"), "DESeq2 ratio method using a prevalence filtered object and a geometric mean that excludes zero counts."),
                                                            tags$li(tags$b("DESeq_poscounts:"), "DESeq2 ratio method with a geometric mean that includes zero counts using the same prevalence filtered object."),
                                                            tags$li(tags$b("library size (RA):"), "Putting all samples to same count level (basically relative abundance)")
                                                    )
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_3')
                                            ),

                                            
                                            wellPanel(
                                                    tags$h4("Calculate size factors"),
                                                    
                                                    textInput(inputId = "prev_SFs", label = "Prevalence for DESeq_gm_exclZero and DESeq_poscounts", value = "60", width = "250px"),
                                                    
                                                    actionButton(inputId = "calcSFs", label = "Calc. size factors")
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Visualize size factor adjustment"),
                                                    
                                                    radioButtons(inputId = "tcaType3", label = "Select type of library size adjustment",
                                                                 choices = c("DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered3", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%"),
                                                    
                                                    actionButton(inputId = "barplotSFs", label = "Visualise library size adjustment"),
                                                    
                                                    tags$h5("NB: Please be patient, calculation of plot can take 20 seconds.")
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            tableOutput(outputId = "overviewView3"),
                                            
                                            #tableOutput(outputId = "explorationViews_2"),
                                           
                                            
                                            plotOutput(outputId = "plotSFs")
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 15 ########################
                   tabPanel(title = "Tax Assignment",
                            sidebarLayout(
                                    sidebarPanel(
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Check taxonomic assignment niveau and whether non-zero counts associate with prevalence.")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_15')
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("The objects to look at here:"),
                                                    
                                                    radioButtons(inputId = "tcaType15", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered15", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%")
                                            ),
                                            
                                            
                                            wellPanel(
                                                    tags$h4("Taxonomic Assignment"),
                                                  
                                                    actionButton(inputId = "calcAssignment", label = "Calc. taxonomic assignment")
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Visualize association of abundance to prevalence"),
                                                    
                                                    actionButton(inputId = "plot_ab_prev", label = "Visualise prevalence to non-zero counts."),
                                                    
                                                    tags$h5("NB: Please be patient, calculation of plot can take a bit.")
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            tableOutput(outputId = "overviewView15"),
                                            
                                            tableOutput(outputId = "table_Assignment"),
                                            
                                            plotOutput(outputId = "plot_Assignment", height = "700px"), #, height = "700px"
                                             
                                            plotOutput(outputId = "plotabPrev", height = "700px")
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 11 ########################
                   tabPanel(title = "Tax_glom",
                            sidebarLayout(
                                    sidebarPanel(
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Option to tax_glom to the chosen taxonomic level")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_11')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    tags$h4("Do tax_glom"),
                                                    
                                                    tags$h5("NB: this action will also refresh all library size adjusted objects if SFs exist!!
                                                            It also resets the filtered taxa."),
                                                    
                                                    radioButtons(inputId = "taxLevel", label = "Select taxonomic level to tax_glom to",
                                                                 choices = c("Species",
                                                                             "Genus",
                                                                             "Family",
                                                                             "Order",
                                                                             "Class",
                                                                             "Phylum"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "NArm", label = "Do you want to remove taxa that are NA at chosen taxonomic level?",
                                                                 choices = c("FALSE",
                                                                             "TRUE"),
                                                                 width = "100%"),
                                                    
                                                    actionButton(inputId = "taxGlom", label = "Perform tax_glom"),
                                                    
                                                    tags$h5("NB: Please be patient, tax_glom can take several minutes")
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            tableOutput(outputId = "overviewView11")
                                            
                                            # tableOutput(outputId = "explorationViews_11")
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 4 ########################
                   tabPanel(title = "Heatmap",
                            sidebarLayout(
                                    sidebarPanel(
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Visualize phyloseq object in a heatmap, and visualize sparsity because zero counts are highlighted by a user-defined color")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_4')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    radioButtons(inputId = "tcaType4", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered4", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%"),
                                                    
                                                    
                                                    textInput(inputId = "max_abundance_for_color", label = "Abundance of max. color", value = "", width = "200px"),
                                                    
                                                    tags$h5("If Abundance of max. color is blank, the max abundance value of the data is used."),
                                                    
                                                    textInput(inputId = "zero_color", label = "Color for zero counts", value = "gray", width = "200px"),
                                                    
                                                    textInput(inputId = "taxa_index_range", label = "taxa index range, e.g. 5:34, leave blank to see all", value = "", width = "200px"),
                                                    
                                                    checkboxInput(inputId = "log4", label = "log abundance", value = FALSE),
                                                    
                                                    actionButton(inputId = "calculateHM", label = "Plot heatmap"),
                                                    
                                                    tags$h5("NB: Please be patient, calculation of plot can take 20 seconds.")
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            # tableOutput(outputId = "explorationViews_3"),
                                            tableOutput(outputId = "overviewView4"),
                                            
                                            plotOutput(outputId = "plotHM")
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 12 ########################
                   tabPanel(title = "Alpha diversity",
                            sidebarLayout(
                                    sidebarPanel(
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Calculate and visualize alpha-diversity.")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_12')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    radioButtons(inputId = "tcaType12", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered12", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%"),
                                                    
                                                    selectInput(inputId = "alpha_measure", label = "Choose the distance measure to use:", choices = c("Observed", "Shannon", "Chao1", "ACE", "Simpson", "InvSimpson", "Fisher"), selected = "Observed"),
                                                    
                                                    checkboxInput(inputId = "rarify", label = "rarify", value = FALSE),
                                                    
     
                                                    
                                                    actionButton(inputId = "calculatealphDiv", label = "calculate alpha diversity"),
                                                    
                                                    tags$h5("NB: Please be patient, calculation can take a bit of time.")
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            tableOutput(outputId = "overviewView12"),
                                            
                                            tableOutput(outputId = "tableAlpha"),
                                            
                                            
                                            plotOutput(outputId = "alphaDiv")
                                            
                                    )

                            )
                   ),
                   ############# Tab 5 ########################
                   tabPanel(title = "Filtering",
                            sidebarLayout(
                                    sidebarPanel(
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Option to prevalence filter the taxa in your phyloseq object for subsequent analyses. Visualize the effect of the filtering on the number of taxa and overall counts in your object. NB: Filtering step has to be done, just choose prevalence = 0 if you do not want to remove any taxa.")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_5')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    
                                                    radioButtons(inputId = "tcaType5", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    textInput(inputId = "filt_prevalence", label = "prevalence", value = "20"),
                                                    
                                                    
                                                    textInput(inputId = "taxa_sums_quantile", label = "taxa_sums quantile in PC", value = "100"),
                                                    
                                                    tags$h5("taxa whose taxa_sums are above this threshold will be kept even if they do not pass prevalence filter"),
                                                    
                                                    wellPanel(
                                                            actionButton(inputId = "prevalenceDistribution", label = "Visualize prevalence distribution"),
                                                            
                                                            actionButton(inputId = "filter", label = "Filter")
                                                    ),
                                                    
                                                    tags$h5("NB: Please be patient, may take ~ 10 seconds.")
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            # tableOutput(outputId = "explorationViews_4"),
                                            tableOutput(outputId = "overviewView5"),
                                            
                                            
                                            plotOutput(outputId = "plotFilter") # height = "1000px"
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 6 ########################
                   tabPanel(title = "Beta diversity",
                            sidebarLayout(
                                    sidebarPanel(
                                            
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Test for overall microbiome differences between the sample groups. 
                                                            A distance (option to choose different distance measures!) is calculated between each two samples in your data. 
                                                            Are samples within the groups closer to each other than between the groups? Do the groups form clusters in an ordination plot?
                                                            (NB: calculations are done on relative abundances of filtered phyloseq object.)")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_6')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    
                                                    radioButtons(inputId = "tcaType6", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered6", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%"),
                                                    
                                                    
                                                    selectInput(inputId = "beta_measure", label = "Choose the distance measure to use:", choices = c("jsd", "bray", "jaccard", "euclidean", "manhattan", "canberra"), selected = "jsd"),

                                                    radioButtons(inputId = "pcoaCorrection", label = "Adjust axes ratio in pcoa?",
                                                                 choices = c("yes",
                                                                             "no"),
                                                                 width = "100%"),
                                                    
                                                    actionButton(inputId = "calcBetaDiversity", label = "Calculate and visualize beta diversity differences")

                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            tableOutput(outputId = "overviewView6"),
                                            
                                            tableOutput(outputId = "adonis"),
                                            
                                            
                                            plotOutput(outputId = "plotPCoA") # height = "1000px"
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 7 ########################
                   tabPanel(title = "Fisher test",
                            sidebarLayout(
                                    sidebarPanel(
                                            
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Find taxa that are more prevalent in one group than in the other.")
                                                    
                                                    ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_7')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    
                                                    tags$h5("NB: the object type only matters for the heatmap here not for the fisher test results."), 
                                                    
                                                    radioButtons(inputId = "tcaType7", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered7", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%"),
                                                    
                                                    actionButton(inputId = "fisher", label = "Calculate prevalence differences using fisher.test"),
                                                    tags$h5(),
                                                    textInput(inputId = "max_abundance_for_colorF", label = "Abundance of max. color", value = "", width = "200px"),
                                                    
                                                    tags$h5("If Abundance of max. color is blank, the max abundance value of the data is used."),
                                                    
                                                    textInput(inputId = "zero_color7", label = "Color for zero counts", value = "gray", width = "200px"),
                                                    
                                                    textInput(inputId = "maxShown7", label = "max number of hits shown in heatmap", value = "40", width = "200px"),
                                                    
                                                    checkboxInput(inputId = "log7", label = "log abundance in heatmap", value = FALSE)
                                                    
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            tableOutput(outputId = "overviewView7"),
  
                                            plotOutput(outputId = "plotFisher", height = "750px"), # height = "1000px"
                                            
                                            tableOutput(outputId = "tableFisher")
                                            
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 8 ########################
                   tabPanel(title = "DESeq2 test",
                            sidebarLayout(
                                    sidebarPanel(
                                            
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Find taxa that are differentially abundant between the groups using DESeq2.")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_8')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    
                                                    radioButtons(inputId = "tcaType8", label = "Select size factor correction",
                                                                 choices = c("DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered8", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%"),
                                                    
                                                    
                                                    actionButton(inputId = "DESeq2", label = "Calculate abundance differences using DESeq2"),
                                                    
                                                    tags$h5(),
                                                    
                                                    textInput(inputId = "max_abundance_for_colorD", label = "Abundance of max. color", value = "", width = "200px"),
                                                    
                                                    tags$h5("If Abundance of max. color is blank, the max abundance value of the data is used."),
                                                    
                                                    textInput(inputId = "zero_color8", label = "Color for zero counts", value = "gray", width = "200px"),
                                                    
                                                    textInput(inputId = "maxShown8", label = "max number of hits shown in heatmap", value = "40", width = "200px"),
                                                    
                                                    checkboxInput(inputId = "log8", label = "log abundance in heatmap", value = FALSE)
                                                    
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            tableOutput(outputId = "overviewView8"),
                                            
                                            plotOutput(outputId = "plotDESeq2", height = "750px"), # height = "1000px"
                                            
                                            tableOutput(outputId = "tableDESeq2")
     
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 9 ########################
                   tabPanel(title = "Wilcoxon test",
                            sidebarLayout(
                                    sidebarPanel(
                                            
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Find taxa that are differentially abundant between the groups using wilcoxon test. NB: Test is done on relative abundances of filtered phyloseq object.")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_9')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    
                                                    radioButtons(inputId = "tcaType9", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered9", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "wilcoxonZeros", label = "Include or exlcude zero counts.",
                                                                 choices = c("include zero counts",
                                                                             "exclude zero counts"),
                                                                 width = "100%"),
                                                    

                                                    actionButton(inputId = "wilcoxon", label = "Calculate abundance differences using wilcoxon"),
                                                    tags$h5(),
                                                    textInput(inputId = "max_abundance_for_colorW", label = "Abundance of max. color", value = "", width = "200px"),
                                                    
                                                    tags$h5("If Abundance of max. color is blank, the max abundance value of the data is used."),
                                                    
                                                    textInput(inputId = "zero_color9", label = "Color for zero counts", value = "gray", width = "200px"),
                                                    
                                                    textInput(inputId = "maxShown9", label = "max number of hits shown in heatmap", value = "40", width = "200px"),
                                                    
                                                    checkboxInput(inputId = "log9", label = "log abundance in heatmap", value = FALSE)
                                                    
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            tableOutput(outputId = "overviewView9"),
                                            
                                            plotOutput(outputId = "plotWilcoxon", height = "750px"), # height = "1000px"
                                            
                                            tableOutput(outputId = "tableWilcoxon")
                                            
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ), 
                   
                   ############# Tab 10 ########################
                   tabPanel(title = "Phylum analysis",
                            sidebarLayout(
                                    sidebarPanel(
                                            
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Compare phylum to phylum ratios. These are independent of compositionality and especially the Firmicutes to Bacteroides ratio has been reported
                                                            extensively in the gut microbiome field in relation to obesity."),
                                                    
                                                    tags$ul(tags$li("Step 1: We calculate the abundance ratios of Firmicutes to all other phyla."),
                                                            tags$li("Step 2: We then calculate a tile plot comparing the ratios of all phyla to each other.")),
                                                    
                                                    tags$h5("NB: in both cases ratios where either the nominator or denominator phylum is absent are set to NA and ignored!
                                                            Wilcoxon test is used for determining whether the ratios are significantly different between the sample groups.")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_10')
                                            ),
                                            
                                            
                                            wellPanel(
                                                    
                                                    tags$h5("Since phylum/phylum analyses are independent of size factor corrections, always the raw counts will be used. But you can decide on whether
                                                            all taxa or only the filtered taxa are included."),
                                                    
                                                    radioButtons(inputId = "tcaType10", label = "Select phyloseq object",
                                                                 choices = c("raw counts"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered10", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%"),
                                                    
                                                    wellPanel(
                                                            
                                                            tags$h5("Step 1: Compare phylum/phylum ratios with boxplots"),
                                                            
                                                            textInput(inputId = "numerator_phylum", label = "Numerator Phylum", value = "Firmicutes", width = "200px"),
                                                            
                                                            textInput(inputId = "denominator_phylum", label = "Denominator Phylum", value = "", width = "200px"),
                                                            
                                                            tags$h5("Leave empty when you want to compare to all other phyla. Otherwise give a comma-separated list without typos."),
                                                            
                                                            actionButton(inputId = "firmicutes", label = "Calculate Firmicutes ratios")
                                                    ),
                                                    
                                                    
                                                    wellPanel(
                                                            tags$h5("Step 2: compare all phylum/phylum ratios"),
                                                            
                                                            actionButton(inputId = "tile", label = "Calculate phylum to phylum tile plot")
                                                    ),
                                                    
                                                    tags$h5("NB: Please be patient, plots can take 10 - 20 seconds.")
                                                    
                                            )
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            tableOutput(outputId = "overviewView10"),
                                            
                                            plotOutput(outputId = "plotPhylum", height = "750px"),
                                            
                                            tableOutput(outputId = "tableFirmicutes")
                                    )
                                    
                            )
                   ),
                   ############# Tab 13 ########################
                   tabPanel(title = "Taxa finder",
                            sidebarLayout(
                                    sidebarPanel(
                                            wellPanel(
                                                    tags$h4("Purpose of this tab:"),
                                                    
                                                    tags$h5("Find specific taxa and compare abundances of those taxa.")
                                                    
                                            ),
                                            
                                            wellPanel(
                                                    tags$h4("Info box"),
                                                    
                                                    textOutput(outputId = 'infoText_13')
                                            ),
                                            
                                            wellPanel(
                                                    
                                                    tags$h4("Select the object to work with on this panel."),
                                                    
                                                    radioButtons(inputId = "tcaType13", label = "Select phyloseq object",
                                                                 choices = c("raw counts",
                                                                             "DESeq_gm_exclZero",
                                                                             "DESeq_poscounts",
                                                                             "library size (RA)"),
                                                                 width = "100%"),
                                                    
                                                    radioButtons(inputId = "filtered13", label = "only filtered taxa?",
                                                                 choices = c("all taxa",
                                                                             "filtered taxa"),
                                                                 width = "100%")
                                                    
                                            ),
                                            
                                            
                                            wellPanel(
                                                    tags$h4("Find taxa"),
                                                    
                                                    radioButtons(inputId = "taxLevelSearch", label = "Select taxonomic level to search in",
                                                                 choices = c("Species",
                                                                             "Genus",
                                                                             "Family",
                                                                             "Order",
                                                                             "Class",
                                                                             "Phylum"),
                                                                 width = "100%"),
                                                    
                                                    textInput(inputId = "searchWord", label = "Search word", value = "Bacteroides"),
                                                    
                                                    actionButton(inputId = "taxaSearch", label = "Perform taxa search")
                                                    
                                                    # tags$h5("NB: Please be patient, tax_glom can take several minutes")
                                            ),
                                            
                                            
                                            wellPanel(
                                                    tags$h4("Plot taxa abundances"),
                                                    
                                                    textInput(inputId = "userTaxa", label = "comma-separated Indexes of Taxa to print", value = ""),
                                                    
                                                    actionButton(inputId = "plotTaxa", label = "Plot abundances of given taxa"),
                                                    
                                                    checkboxInput(inputId = "log", label = "log abundance", value = FALSE),
                                                    
                                                    checkboxInput(inputId = "pool", label = "pool abundances", value = FALSE)
                                                    
                                                    # tags$h5("NB: Please be patient, tax_glom can take several minutes")
                                            )
                                            
                                            
                                            
                                    ),
                                    mainPanel(
                                            # tags$h2("The current phyloseq object"),
                                            
                                            tableOutput(outputId = "overviewView13"),
                                            
                                            tableOutput(outputId = "foundTaxa"),
                                            
                                            tableOutput(outputId = "pValsTaxa"),
                                            
                                            plotOutput(outputId = "plotTaxaAbundances")
                                            
                                    )
                                    # textOutput(outputId = 'explorationViews'))
                            )
                   ),
                   
                   ############# Tab 14 ########################
                   {
                           tabPanel(title = "Taxa ratio",
                                    sidebarLayout(
                                            sidebarPanel(
                                                    wellPanel(
                                                            tags$h4("Purpose of this tab:"),
                                                            
                                                            tags$h5("Compare the ratios of two user defined taxa.")
                                                            
                                                    ),
                                                    
                                                    wellPanel(
                                                            tags$h4("Info box"),
                                                            
                                                            textOutput(outputId = 'infoText_14')
                                                    ),
                                                    
                                                    
                                                    wellPanel(
                                                            tags$h4("plot taxa Abundance ratios"),
                                                            
                                                            radioButtons(inputId = "tcaType14", label = "Select phyloseq object",
                                                                         choices = c("raw counts",
                                                                                     "DESeq_gm_exclZero",
                                                                                     "DESeq_poscounts",
                                                                                     "library size (RA)"),
                                                                         width = "100%"),
                                                            
                                                            radioButtons(inputId = "filtered14", label = "only filtered taxa?",
                                                                         choices = c("all taxa",
                                                                                     "filtered taxa"),
                                                                         width = "100%"),
                                                            
                                                            
                                                            textInput(inputId = "numerator", label = "Index of numerator", value = ""),
                                                            
                                                            textInput(inputId = "denominator", label = "Index of denominator", value = ""),
                                                            
                                                            checkboxInput(inputId = "logRatio", label = "log Ratios", value = FALSE),
                                                            
                                                            actionButton(inputId = "plotRatioTaxa", label = "plot abundance ratios of the taxa")
                                                            
                                                            # tags$h5("NB: Please be patient, tax_glom can take several minutes")
                                                    )
                                                    
                                                    
                                                    
                                            ),
                                            mainPanel(
                                                    # tags$h2("The current phyloseq object"),
                                                    
                                                    tableOutput(outputId = "overviewView14"),
                                                    
                                                    tableOutput(outputId = "pValsRatioTaxa"),
                                                    
                                                    plotOutput(outputId = "plotTaxaRatios")
                                                    
                                            )
                                            # textOutput(outputId = 'explorationViews'))
                                    )
                           )
                   }
            
        )
)







server <- function(input, output, session){
        
        
        
        ################################ Overalls ############################################
        # - generate reactive Values for entire app -
        
        
        
        rv <- reactiveValues(viewItem = NULL, overviewView = NULL, ps = NULL, infoText = NULL, infoText_2 = NULL, phyla = NULL, phylum_colors = NULL, Tr = NULL, Tr_SFs = NULL,
                             ps_tca = NULL, ps_tca_poscounts = NULL, SFs = NULL, SFs_poscounts = NULL, SFs_RA = NULL, ps_RA = NULL, Tr_HM = NULL, 
                             filteredKeepTaxa = NULL, tableAssignment = NULL, Tr_Assignment = NULL, Tr_ab_prev = NULL, TrL_filter = NULL, Tr_PCoA = NULL, adonis = NULL,
                             Tr_fisher = NULL, fisherTable = NULL, Tr_DESeq2 = NULL, DESeq2Table = NULL, Tr_wilcoxon = NULL, wilcoxonTable = NULL,
                             firmicutesTable = NULL, Tr_phylum = NULL, Tr_alpha = NULL, alphaTable = NULL, group_var = NULL, group_var_levels = NULL,
                             color_levels = NULL, taxaFindTable = NULL, infoText_13 = NULL, Tr_taxaAbundance = NULL,
                             tablepValsTaxa = NULL, Tr_taxaRatios = NULL, tablepValsTaxaRatios = NULL,
                             symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

        
        ################################ Tab 1 ############################################
        # - Tab1: Parameters -
        observe({
                rv$group_var <- input$group_var
                grp1 <- input$grp1
                grp2 <- input$grp2
                group_var_levels <- c(grp1, grp2)
                rv$group_var_levels <- group_var_levels
                grp1_color <- input$grp1_color
                grp2_color <- input$grp2_color
                color_levels <- c(grp1_color, grp2_color)
                names(color_levels) <- rv$group_var_levels
                rv$color_levels <- color_levels
        })
        
        
        observe({
                ps <- rv$ps
                ps_tca <- rv$ps_tca 
                ps_tca_poscounts <- rv$ps_tca_poscounts
                ps_RA <- rv$ps_RA
                filteredTaxa <- rv$filteredKeepTaxa
                if (is.null(filteredTaxa)){
                        filteredTaxa <- "NULL"
                } else {
                        filteredTaxa <- paste0(length(filteredTaxa), " taxa")
                }
                
                Items <- c("ps: raw_counts",
                           "ps: DESeq_gm_exclZero",
                           "ps: DESeq_poscounts",
                           "ps: relative Abundance",
                           "filtered taxa")
                Status <- c(output_ps(ps), output_ps(ps_tca), output_ps(ps_tca_poscounts),
                            output_ps(ps_RA), filteredTaxa)
                
                rv$overviewView <- data.frame(Item = Items, Status = Status)
                
        })
        
        
        # -- render the overviewViews table --
        output$overviewView1 <- renderTable({
                rv$overviewView
        }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
        # ----
        # --
        
        
        
        
        ################################ Tab 2 ############################################
        {    # NB: I add these {} so I can make the entire Tab small: makes it easier to work on individual tabs    
                # - Tab2: Load and explore -
                # -- output infoText --
                output$infoText_2 <- renderText({
                        rv$infoText_2
                })
                # ----
        
                # -- render the overviewViews table --
                output$overviewView2 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                # -- Load phyloseq object --
                observeEvent(input$loadPhyseq, {
                        
                        rv$ps <- NULL
                        rv$ps_tca <- NULL
                        rv$ps_tca_poscounts <- NULL
                        rv$ps_RA <- NULL
                        rv$Tr <- NULL
                        rv$infoText_2 <- NULL
                        rv$phyla <- NULL
                        rv$phylum_colors <- NULL
                        rv$filteredKeepTaxa <- NULL
        
                        inFile <- input$loadPhyseq
                        
                        if(is.null(inFile)){
                                rv$infoText_2 <- "No correct rds file was selected. Please try again."
                                return()
                        }
                        
                        if (!grepl(pattern = ".rds$", inFile$datapath)) {
                                rv$infoText_2 <- "Sorry. The chosen file was not a .rds file."
                                return()
                                
                        }
                        
                        
                        ps <- readRDS(file = inFile$datapath)
                        
                        
                        if (class(ps) != "phyloseq"){
                                rv$infoText_2 <- "The chosen .rds file did not contain a phyloseq object. Please choose phyloseq_object_Zeller_CRC for this tutorial."
                                return()
                                
                        }
                        
                        
                        
                        # - added to remove taxa that are not present in a single sample directly -
                        keepTaxa <- taxa_names(ps)[taxa_sums(ps) > 0]
                        ps <- phyloseq::prune_taxa(keepTaxa, ps)
                        # --
                        
                        
                        
                        rv$ps <- ps
                        rv$infoText_2 <- paste("Uploaded phyloseq object with ", nsamples(ps), " samples and ", ntaxa(ps), " taxa.", sep = "")
        
                })
                # ----
        
                
        
                # # -- render the explorationViews table --
                # output$explorationViews <- renderTable({
                #         if(!is.null(rv$ps)){
                #                 psShow <- capture.output(rv$ps)
                #                 psShow <- as.data.frame(psShow, nrow = 4)
                #                 psShow <- psShow[2:4, , drop = FALSE]
                #                 colnames(psShow) <- "object component"
                #                 psShow <- tidyr::separate(psShow, col = "object component", into = c("object component", "dimension"), sep = ":")
                #                 psShow
                #         } else {
                #                 NULL
                #         }
                # }, sanitize.text.function = function(x) x, caption = "Loaded phyloseq object", caption.placement = getOption("xtable.caption.placement", "top"))
                # # ----
        
        
        
                # -- visualise physeq components --
                observeEvent(input$checkPhyseq, {
                        
                        rv$phyla <- NULL
                        rv$infoText_2 <- NULL
                        rv$Tr <- NULL
                        
                        if (is.null(rv$ps)) {
                                rv$infoText_2 <- "Sorry. There is no phyloseq object loaded."
                                return()
                        }
                        
                        # -- evaluate sample and taxaRanges --
                        sampleRange <- input$sampleRange
                        if (is.character(sampleRange)){
                                sampleRange <- strsplit(x = input$sampleRange, split = ":")
                        } else {
                                sampleRange <- list(paste0(1, ":", nsamples(rv$ps)))
                        }
                        sampleRange <- strsplit(x = input$sampleRange, split = ":")
                        firstSample <- as.numeric(sapply(sampleRange, `[`, 1))
                        lastSample <- as.numeric(sapply(sampleRange, `[`, 2))
                        taxaRange <- strsplit(x = input$taxaRange, split = ":")
                        firstTaxon <- as.numeric(sapply(taxaRange, `[`, 1))
                        lastTaxon <- as.numeric(sapply(taxaRange, `[`, 2))
                        
                        if (is.na(firstSample) || is.na(lastSample) || !is.numeric(firstSample) || !is.numeric(lastSample) ||
                            firstSample < 1 || lastSample > nsamples(rv$ps) || firstSample > lastSample) {
                                
                                firstSample <- 1
                                lastSample <- nsamples(rv$ps)
                        }
                        
                        if (is.na(firstTaxon) || is.na(lastTaxon) || !is.numeric(firstTaxon) || !is.numeric(lastTaxon) ||
                            firstTaxon < 1 || lastTaxon > ntaxa(rv$ps) || firstTaxon > lastTaxon) {
                                
                                firstTaxon <- 1
                                lastTaxon <- ntaxa(rv$ps)
                        }
                        # ----
                        
                        
                        if (input$component == "sample_data"){
                                
                                SD <- as(sample_data(rv$ps), "data.frame")
                                
                                SD <- SD[firstSample:lastSample, ]
                                
                                rv$phyla <- SD
                                rv$infoText_2 <- "Sample data has been extracted and will be shown."
                                
                        } else if (input$component == "otu_table"){
                                
                                if (taxa_are_rows(rv$ps)){
                                        OTU <- round(as(otu_table(rv$ps), "matrix"))
                                } else {
                                        OTU <- round(t(as(otu_table(rv$ps), "matrix")))
                                }
                                
                                OTU <- as.data.frame(OTU)
                                
                                OTU <- OTU[firstTaxon:lastTaxon, firstSample:lastSample]
                                
                                
                                rv$phyla <- OTU
                                rv$infoText_2 <- "count table has been extracted and will be shown."
          
                        } else if (input$component == "tax_table"){
                                TT <- as.data.frame(unclass(tax_table(rv$ps)))
                                
                                TT <- TT[firstTaxon:lastTaxon,]
                                rv$phyla <- TT
                                rv$infoText_2 <- "taxa table has been extracted and will be shown."
                                
                        }
                        
                         
                })
                # ----
        

                
                # -- calc phyla distribution and colors --
                observeEvent(input$calcPhyla, {
                        
                        rv$phylum_colors <- NULL
                        rv$phyla <- NULL
                        rv$infoText_2 <- NULL
                        rv$Tr <- NULL
                        
                        if (is.null(rv$ps)) {
                                rv$infoText_2 <- "Sorry. There is no phyloseq object loaded."
                                return()
                        }
                        
                        
                        if (input$tcaType2 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType2 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType2 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType2 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_2 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_2 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        if (input$filtered2 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_2 <- "Sorry. You have to do the filtering first."
                                        return()
                                        
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                                
                        }
                        
                        Phyla <- check_phyla_distribution_NA(ps)
                        PhylaForColor <- check_phyla_distribution(ps)
                        
                        if (nrow(PhylaForColor) < 16) {
                                phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), QuantColors15)
                        } else {
                                phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), c(QuantColors15, viridis(nrow(PhylaForColor)-15)))
                        }
                        
                        rv$phylum_colors <- phylum_colors
                        Phyla <- dplyr::select(Phyla, Phylum:PC_of_counts, mean_prevalence_in_PC)
                        rv$phyla <- Phyla
                        rv$infoText_2 <- "Phylum distribution has been calculated and phylum_colors have been determined." 
                })
                # ----
        
                
                
                # -- render the phylaViews table --
                output$phylaViews <- renderTable({
                        if(!is.null(rv$phyla)){
                                phylaShow <- rv$phyla
                                phylaShow
                        } else {
                                NULL
                        }
                }, rownames = TRUE, digits = 0, sanitize.text.function = function(x) x, caption = "Phyla distribution or phyloseq component", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        

                # -- generate a barplot of ps --
                observeEvent(input$barplotPS, {
                        
                        rv$Tr <- NULL
                        rv$infoText_2 <- NULL
                        rv$phyla <- NULL
                        
                        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                                rv$infoText_2 <- "Sorry. Missing phyloseq object or phylum_colors. Please do these steps first."
                                return()
                                
                        }
                        
                        phylum_colors <- rv$phylum_colors
                        
                        if (input$tcaType2 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType2 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType2 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType2 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_2 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_2 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        if (input$filtered2 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_2 <- "Sorry. You have to do the filtering first."
                                        return()
                                        
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                                
                        }
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_2 <- error_message
                                return()
                        }
                        # ------
                        
                        
                        psP <- phyloseq::tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
                        
                        Tr <- plot_sample_bars(physeq = psP, x = "Sample", y = "Abundance", group_var, color_levels, fill = "Phylum",
                                               color_sample_names = TRUE, col_vec = phylum_colors, facet_grid = NULL, order_by_firmicutes = FALSE)
 
                        
                        rv$Tr <- Tr
                        rv$infoText_2 <- "Barplot of samples has been generated." 
                })
                # ----
        
        
                # -- output plot samplesPS --
                output$samplesPS <- renderPlot({
                        if (is.null(rv$Tr)) {
                                rv$Tr
                        } else {
                                rv$Tr
                        }
                }, 
                height = 610)
                # ----
                # --
        }

        ################################ Tab 3 ############################################
        {
                # - Tab3: library size adjustment -
                # -- output infoText --
                output$infoText_3 <- renderText({
                        rv$infoText_3
                })
                # ----
        
                
                # -- render the overviewViews table --
                output$overviewView3 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                # # -- render the explorationViews table --
                # output$explorationViews_2 <- renderTable({
                #         if(!is.null(rv$ps)){
                #                 psShow <- capture.output(rv$ps)
                #                 psShow <- as.data.frame(psShow, nrow = 4)
                #                 psShow <- psShow[2:4, , drop = FALSE]
                #                 colnames(psShow) <- "object component"
                #                 psShow <- tidyr::separate(psShow, col = "object component", into = c("object component", "dimension"), sep = ":")
                #                 psShow
                #         } else {
                #                 NULL
                #         }
                # }, sanitize.text.function = function(x) x, caption = "Loaded phyloseq object", caption.placement = getOption("xtable.caption.placement", "top"))
                # # ----
        
        
        
                # -- calc size factors --
                observeEvent(input$calcSFs, {
                        
                        rv$infoText_3 <- NULL
                        rv$ps_tca <- NULL
                        rv$ps_tca_poscounts <- NULL
                        rv$ps_RA <- NULL
                        rv$SFs <- NULL
                        rv$SFs_poscounts <- NULL
                        rv$SFs_RA <- NULL
                        
                        
                        if (is.null(rv$ps)) {
                                rv$infoText_3 <- "Sorry. You need to load a phyloseq object first (Tab 2)."
                                return()
                        }
                        
                        prevalence_for_sf<- as.numeric(input$prev_SFs)
                        
                        if (is.na(prevalence_for_sf) || prevalence_for_sf > 100 || prevalence_for_sf < 0){
                                rv$infoText_3 <- "Sorry. Prevalence for DESeq_gm_exclZero must be between 0 and 100. Please change your input."
                                return()
                                
                        }
                        
                        min_obs <- 0L
                        
                        ps_sf_filt <- phyloseq::filter_taxa(rv$ps, function(x){(sum(x > min_obs) > (prevalence_for_sf/100)*length(x))}, prune = TRUE)
                        
                        SFs <- calc_SFs(physeq = ps_sf_filt)
                        
                        
                        group_var <- rv$group_var
                        
                        if(! group_var %in% colnames(sample_data(rv$ps))) {
                                rv$infoText_3 <- "The given group_var is not a variable in the sample data of the loaded phyloseq object. A correct group variable is needed for size factor calculation using poscount method"
                                return()
                        }
                        
                        SFs_poscounts <- calc_SFs_DESeq(ps_sf_filt, type = "poscounts", group_var = group_var)
                        
                        SFs_RA <- sample_sums(rv$ps)/gm_own(sample_sums(rv$ps), zeros.count = FALSE)
                        # If you want real relative abundances instead use:
                        # SFs_RA <- sample_sums(rv$ps) 
                        # But NB: will cause trouble with alpha diversity and also DESeq2!
                        
                        rv$SFs <- SFs
                        rv$SFs_poscounts <- SFs_poscounts
                        rv$SFs_RA <- SFs_RA
                        
                        library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs) 
                        rv$ps_tca <- library_size_adjust_list[[1]]
                        
                        library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs_poscounts) 
                        rv$ps_tca_poscounts <- library_size_adjust_list[[1]]
                        
                        library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = SFs_RA) 
                        rv$ps_RA <- library_size_adjust_list[[1]]
                        
                        rv$infoText_3 <- "Library size adjusted phyloseq objects have been calculated." 
                })
                # ----
        

                
                # -- generate a barplot before after of selected library size adjustment --
                observeEvent(input$barplotSFs, {
                        
                        rv$infoText_3 <- NULL
                        rv$Tr_SFs <- NULL
                        
                        if (is.null(rv$ps)) {
                                rv$infoText_3 <- "Sorry. No phyloseq object is loaded."
                                return()
                        }
                        
                        ps <- rv$ps
                        
                        if (is.null(rv$phylum_colors)) {
                                rv$infoText_3 <- "Sorry. phylum_colors have not yet been calculated."
                                return()
                                
                        }
                        
                        phylum_colors <- rv$phylum_colors
                        
                        
                        if (input$tcaType3 == "raw counts"){
                                plot_ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType3 == "DESeq_gm_exclZero"){
                                plot_ps <- rv$ps_tca
                        } else if (input$tcaType3 == "DESeq_poscounts") {
                                plot_ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType3 == "library size (RA)") {
                                plot_ps <- rv$ps_RA
                        } else {
                                rv$infoText_3 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(plot_ps)) {
                                rv$infoText_3 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        if (input$filtered3 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_3 <- "Sorry. You have to do the filtering first."
                                        return()
                                        
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                                ps_plot <- prune_taxa(keepTaxa, ps)
                                
                        }
                        
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_3 <- error_message
                                return()
                        }
                        # ------
                        

                        
                        psP <- phyloseq::tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
                        ps_tcaP <- phyloseq::tax_glom(plot_ps, taxrank = "Phylum", NArm = FALSE)
                        rv$Tr_SFs <- plot_sample_bars_compare(physeq = psP, physeq2 = ps_tcaP, x = "Sample", y = "Abundance", group_var = group_var, color_levels = color_levels, color_sample_names = TRUE, fill = "Phylum", col_vec = phylum_colors, order_by_raw_counts = TRUE)
                        
                        
                        rv$infoText_3 <- "Plot illustrating library size adjustment hast been generated." 
                })
                # ----
        
                
        
                # -- output plot samplesPS --
                output$plotSFs <- renderPlot({
                        if (is.null(rv$Tr_SFs)) {
                                rv$Tr_SFs
                        } else {
                                rv$Tr_SFs
                        }
                },
                height = 700)
                # ----
                # --
        
        }
        
        ################################ Tab 15 ############################################
        {
                # - Tab15: taxonomic assignment -
                # -- output infoText --
                output$infoText_15 <- renderText({
                        rv$infoText_15
                })
                # ----
                
                
                # -- render the overviewViews table --
                output$overviewView15 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                # -- calc taxonomic assignment --
                observeEvent(input$calcAssignment, {
                        
                        rv$infoText_15 <- NULL
                        rv$tableAssignment <- NULL
                        rv$Tr_Assignment <- NULL
                        rv$Tr_ab_prev <- NULL
                        
                        if (input$tcaType15 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType15 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType15 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType15 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_15 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_15 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        if (input$filtered15 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_15 <- "Sorry. You have to do the filtering first."
                                        return()
                                        
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                                
                        }
                        
                        
                        assignment_distribution <- get_assignemnt_distribution(ps)
                        
                        assignment_distribution <- cbind(Level = rownames(assignment_distribution), assignment_distribution)
                        
                        rv$tableAssignment <- assignment_distribution
                        
                        assign_vs_ab <- check_assignment_vs_abundance(ps)
                        assign_vs_prev <- check_assignment_vs_prevalence(ps)
                        
                        rv$Tr_Assignment <- list(assign_vs_prev[[2]], assign_vs_ab[[2]])
                        
                        rv$infoText_15 <- "Taxonomic assignment niveaus have been calculated."
                        
                })
                # ----
                
                
                # -- render the assignment table --
                output$table_Assignment <- renderTable({
                        if(!is.null(rv$tableAssignment)){
                                rv$tableAssignment
                        } else {
                                NULL
                        }
                }, digits = 2, sanitize.text.function = function(x) x, caption = "Taxonomic assignment levels:", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                # -- output plot of Assignment curves --
                output$plot_Assignment <- renderPlot({
                        if (is.null(rv$Tr_Assignment)) {
                                rv$Tr_Assignment
                        } else {
                                do.call("grid.arrange", c(rv$Tr_Assignment, nrow = 2))
                        }
                })
                # ----
                
                
                # -- plot prevalence to non-zero abundance --
                observeEvent(input$plot_ab_prev, {
                        
                        rv$infoText_15 <- NULL
                        rv$tableAssignment <- NULL
                        rv$Tr_Assignment <- NULL
                        rv$Tr_ab_prev <- NULL
                        
                        if (input$tcaType15 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType15 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType15 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType15 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_15 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_15 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        if (input$filtered15 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_15 <- "Sorry. You have to do the filtering first."
                                        return()
                                        
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                                
                        }
                        
                        
                        if (is.null(rv$phylum_colors)) {
                                rv$infoText_15 <- "Sorry. You need to calculate phylum_colors first."
                                return()
                                
                        }
                        
                        phylum_colors <- rv$phylum_colors
                        
                        TrrList <- plot_correlations_abundance_prev_sparsity(physeq = ps, col = "Phylum", col_vec = phylum_colors)
                        
                        rv$Tr_ab_prev <- list(TrrList[[3]], TrrList[[4]])
                        
                        rv$infoText_15 <- "Plots of prevalence to non-zero abundances have been calculated."
                        
                })
                # ----
                
                
                
                # -- output plot of Assignment curves --
                output$plotabPrev <- renderPlot({
                        if (is.null(rv$Tr_ab_prev)) {
                                rv$Tr_ab_prev
                        } else {
                                do.call("grid.arrange", c(rv$Tr_ab_prev, ncol = 2))
                        }
                })
                # ----
                
                # --
                
        }
        
        ################################ Tab 11 ############################################
        {
                # - Tab11: Tax_glom -
                # -- output infoText --
                output$infoText_11 <- renderText({
                        rv$infoText_11
                })
                # ----
        
        
                # -- render the overviewViews table --
                output$overviewView11 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
                # -- tax_glom --
                observeEvent(input$taxGlom, {
                        
                        rv$infoText_11 <- NULL
                        rv$filteredKeepTaxa <- NULL
                        
                        if (is.null(rv$ps)) {
                                rv$infoText_11 <- "Sorry. No phyloseq object is loaded."
                                return()
                        }
                        
                        ps <- rv$ps
                        
                        tax_level <- input$taxLevel
                        NArmal <- as.logical(input$NArm)
                        
                        ps <- phyloseq::tax_glom(ps, taxrank = tax_level, NArm = NArmal)
                        
                        rv$ps <- ps
                        
                        # --- adjust size factor corrected objects if SFs have already been calculated ---
                        
                        if (!is.null(rv$SFs) && !is.null(rv$SFs_poscounts) && !is.null(rv$SFs_RA)) {
                                
                                library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = rv$SFs) 
                                rv$ps_tca <- library_size_adjust_list[[1]]
                                
                                library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = rv$SFs_poscounts) 
                                rv$ps_tca_poscounts <- library_size_adjust_list[[1]]
                                
                                library_size_adjust_list <- simply_adjust_LS(physeq = rv$ps, SFs = rv$SFs_RA) 
                                rv$ps_RA <- library_size_adjust_list[[1]]
                                
                                rv$infoText_11 <- paste("Tax_glom has been performed on original object and size factor adjusted objects!", sep = "")
                                
                        } else {
                                
                                rv$infoText_11 <- paste("Tax_glom has been performed on original object. No size factors had been calculated yet.", sep = "")
                                
                        }
                        
                        # ------
        
                })
                # ----
        }
        ################################ Tab 4 ############################################
        {
                # - Tab4: heatmap/sparsity -
                # -- output infoText --
                output$infoText_4 <- renderText({
                        rv$infoText_4
                })
                # ----
                
                
                # -- render the overviewViews table --
                output$overviewView4 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
                
                # -- calc heatmap --
                observeEvent(input$calculateHM, {
                        
                        rv$infoText_4 <- NULL
                        rv$Tr_HM <- NULL
                        
                        
                        if (input$tcaType4 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType4 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType4 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType4 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_4 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_4 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        if (input$filtered4 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_4 <- "Sorry. You have to do the filtering first."
                                        return()
                                        
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                                
                        }
                        
                        zero_color <- input$zero_color
                        
                        if (!areColors(zero_color)){
                                rv$infoText_4 <- "Sorry. The color for the zero counts is not an R color, please change."
                                return()
                                
                        }
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_4 <- error_message
                                return()
                        }
                        # ------
                        
                        
                        # --- taxa_range option ---
                        taxa_index_range <- input$taxa_index_range
                        
                        if (taxa_index_range != ""){
                                
                                taxa_index_range <- strsplit(x = taxa_index_range, split = ":")
                                taxa_index_range <- round(as.numeric(unlist(taxa_index_range))[1:2])
                                if (any(is.na(taxa_index_range))){
                                        rv$infoText_4 <- "Sorry couldn't read taxa_index_range, please change input."
                                        return()
                                }
                                start_index <- taxa_index_range[1]
                                end_index <- taxa_index_range[2]
                                
                                if (start_index > end_index || start_index < 1 || end_index > ntaxa(ps)){
                                        rv$infoText_4 <- "Sorry taxa_index_range did not fit to ps or start index was bigger than end index, please change input."
                                        return()
                                }
                                
                        } else {
                                start_index <- 1
                                end_index <- ntaxa(ps)
                        }
                        
                        keepTaxa <- taxa_names(ps)[start_index:end_index]
                        ps <- prune_taxa(keepTaxa, ps)
                        # ------
                        
                        
                        
                        df <- as.data.frame(as(tax_table(ps), "matrix"))
                        taxa_annotation <- get_taxon_names_plusTL(df)
                        taxa_annotation <- strsplit(taxa_annotation, split = "/")
                        taxa_annotation <- sapply(taxa_annotation, `[`, 1)
                        taxa_annotation <- make.unique(taxa_annotation)
                        
                        sample_colors <- list(color_levels)
                        names(sample_colors) <- group_var
                        
                        max_abundance_for_color <- input$max_abundance_for_color
                        
                        if (max_abundance_for_color == ""){
                                max_abundance_for_color <- max(otu_table(ps))
                        } else {
                                max_abundance_for_color <- as.numeric(max_abundance_for_color)
                                
                        }
                        
                        if (is.na(max_abundance_for_color) || max_abundance_for_color < 0 || max_abundance_for_color > max(otu_table(ps), na.rm = TRUE)) {
                                rv$infoText_4 <- "Sorry. max_abundance_for_color either not numeric, too small or too high:). Change it."
                                return()
                                
                        }
                                
                        if (input$log4){
                                logger <- TRUE 
                        } else {
                                logger <- FALSE
                        }
                        
                        rv$Tr_HM <- plot_heatmap_physeq(physeq = ps, sample_colors = sample_colors, taxa_info_df = NULL, taxa_colors = NULL,
                                                 taxa_annotation = taxa_annotation, max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1),
                                                 zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = logger, drop_color_levels = TRUE,
                                                 border_color = NA,
                                                 cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE,
                                                 annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 16,
                                                 fontsize_row = 8, fontsize_col = 8, fontsize_number = 12)
        
                        
                        sparsity <- round(100*sum(otu_table(ps) == 0)/length(otu_table(ps)), 2)
                        
                        rv$infoText_4 <- paste("Heatmap has been generated. There are ", ntaxa(ps), " taxa in the plotted data. Overall sparsity of the count table is ", sparsity, "%.", sep = "") 
                })
                # ----

                
        
                # -- output plot samplesPS --
                output$plotHM <- renderPlot({
                        if (is.null(rv$Tr_HM)) {
                                NULL
                        } else {
                                grid::grid.newpage()
                                grid::grid.draw(rv$Tr_HM)
                        }
                },
                height = 1200)
                # ----
                # --
        }
        

        ################################ Tab 12 ############################################
        {
                # - Tab12: alpha diversity -
                # -- output infoText --
                output$infoText_12 <- renderText({
                        rv$infoText_12
                })
                # ----
                
                
                # -- render the overviewViews table --
                output$overviewView12 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                
                # -- calc alpha Diversity --
                observeEvent(input$calculatealphDiv, {
                        
                        rv$infoText_12 <- NULL
                        
                        rv$Tr_alpha <- NULL
                        
                        rv$alphaTable <- NULL
                        
                        
                        if (input$tcaType12 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType12 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType12 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType12 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_12 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_12 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        if (input$filtered12 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_12 <- "Sorry. You have to do the filtering first."
                                        return()
                                        
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                                
                        }

                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_12 <- error_message
                                return()
                        }
                        # ------
                        
                        
                        
                        # --- NB: estimateR.default, in phyloseq::estimate_richness accepts only integers therefore: ---
                        OTU <- as(otu_table(ps), "matrix")
                        if (mode(OTU) != "integer"){
                                OTU <- round(OTU, 1)
                                mode(OTU) <- "integer"
                                otu_table(ps) <- otu_table(OTU, taxa_are_rows = taxa_are_rows(ps))
                        }
                        # -------
                        
                        
                        if (input$rarify) {
                                rare_level <- min(sample_sums(ps))
                                count_table_rare <- vegan::rrarefy(as(otu_table(ps), "matrix"), sample = rare_level)
                                
                                otu_table(ps) <- otu_table(count_table_rare, taxa_are_rows = taxa_are_rows(ps))
                        }
                        
                        
                        alpha_div_measures <- input$alpha_measure
                        
                        DF_alpha_list <- calc_alphadiv_plusLmResids(physeq = ps, measures = alpha_div_measures, group_var = group_var,
                                                                    compare = names(color_levels))
                        
                        lm_fitlist <- DF_alpha_list[[2]]
                        DF_alpha <- DF_alpha_list[[1]]
                        
                        alpha_div_pVals <- calc_pVals_alphdiv(DF_alpha = DF_alpha, measures = alpha_div_measures, group_var = group_var, compare = names(color_levels), test = "t.test")
                        

                        alpha_div_boxplots <- boxplots_alphdiv(DF_alpha = DF_alpha, measures = alpha_div_measures, group_var = group_var, shape = NULL, color_levels = color_levels, test = "t.test", hide.ns = FALSE)

                        alpha_div_lmPlots <- lmPlots_alphdiv(DF_alpha = DF_alpha, lm_fitlist = lm_fitlist, measures = alpha_div_measures, group_var = group_var, shape = NULL, color_levels = color_levels, test = "t.test")

                        
                        TrList <- c(alpha_div_boxplots, alpha_div_lmPlots)
                        TrList <- TrList[order(names(TrList))] 
                        
                        
                        rv$Tr_alpha <- TrList
                        
                        rv$alphaTable <- alpha_div_pVals
                        
                        
                        rv$infoText_12 <- "alpha diversity analysis has been performed."
                        
                })
                # ----
                
                
                
                # -- render the alpha table --
                output$tableAlpha <- renderTable({
                        if(!is.null(rv$alphaTable)){
                                alphaShow <- rv$alphaTable
                                alphaShow
                        } else {
                                NULL
                        }
                }, digits = 4, sanitize.text.function = function(x) x, caption = "alpha diversity test", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                
                # -- output alpha diversity plot --
                output$alphaDiv <- renderPlot({
                        if (is.null(rv$Tr_alpha)) {
                                NULL
                        } else {
                                do.call("grid.arrange", c(rv$Tr_alpha, ncol = 3))
                        }
                },
                height = 500)
                # ----
                
        }
        ################################ Tab 5 ############################################
        {
                # - Tab5: filtering -
                # -- output infoText --
                output$infoText_5 <- renderText({
                        rv$infoText_5
                })
                # ----
        
                
        
                # -- render the overviewViews table --
                output$overviewView5 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----

        
        
                # -- show prevalence distribution --
                observeEvent(input$prevalenceDistribution, {
                        
                        rv$infoText_5 <- NULL
                        
                        rv$TrL_filter <- NULL
                        
                        if (input$tcaType5 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType5 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType5 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType5 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_5 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_5 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        prevalence <- as.numeric(input$filt_prevalence)
                        
                        if (is.na(prevalence) || prevalence > 100 || prevalence < 0){
                                rv$infoText_3 <- "Sorry. Prevalence for filtering must be between 0 and 100. Please change your input."
                                return()
                                
                        }
                        
                        TrList <- plot_ab_pev_distributions(ps, prevalence = prevalence)
                        
                        rv$TrL_filter <- list(TrList[[2]], TrList[[4]])
                        
                        rv$infoText_5 <- "Filter test has been performed."
                        
                })
                # ----
        
                
        
                # -- do the filtering --
                observeEvent(input$filter, {
                        
                        rv$infoText_5 <- NULL
                        rv$TrL_filter <- NULL
                        rv$filteredKeepTaxa <- NULL
                        
                        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                                rv$infoText_5 <- "Sorry. You need a loaded phyloseq object and phylum_colors defined."
                                return()
                                
                        }
                        
                        if (input$tcaType5 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType5 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType5 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType5 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_5 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_5 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        prevalence <- as.numeric(input$filt_prevalence)
                        
                        if (is.na(prevalence) || prevalence > 100 || prevalence < 0){
                                rv$infoText_5 <- "Sorry. Prevalence for filtering must be between 0 and 100. Please change your input."
                                return()
                                
                        }
                        
                        taxa_sums_quantile <- as.numeric(input$taxa_sums_quantile)
                        
                        if (is.na(taxa_sums_quantile) || taxa_sums_quantile > 100 || taxa_sums_quantile < 0){
                                rv$infoText_5 <- "Sorry. Prevalence for taxa_sums_quantile must be between 0 and 100. Please change your input."
                                return()
                                
                        }
                        
                        
                        min_obs <- 0L
                        
                        ps_Filt <- phyloseq::filter_taxa(ps, function(x){
                                (sum(x > min_obs) > (prevalence/100)*length(x)) || 
                                        (sum(x) > quantile(taxa_sums(ps), probs = taxa_sums_quantile/100))
                        }, prune = TRUE)
                        
                        keepTaxa <- taxa_names(ps_Filt)
                        
                        rv$filteredKeepTaxa <- keepTaxa
                        
                        
                        rv$TrL_filter <- visualize_filtering(physeq = ps, prevalence = prevalence, taxa_sums_quantile = taxa_sums_quantile, phylum_colors = rv$phylum_colors)
                
                        
                        rv$infoText_5 <- "ps_filt has been generated. Filtering done." 
                })
                # ----
        
        
        
                # -- output plot filter --
                output$plotFilter <- renderPlot({
                        if (is.null(rv$TrL_filter)) {
                                rv$TrL_filter
                        } else {
                                grid.arrange(rv$TrL_filter[[1]], rv$TrL_filter[[2]], nrow = 2)
                        }
                },
                height = 800)
                # --

        }                
        
        ################################ Tab 6 ############################################
        {
                # - Tab6: beta diversity -
                # -- output infoText --
                output$infoText_6 <- renderText({
                        rv$infoText_6
                })
                # ----
                
                # -- render the overviewViews table --
                output$overviewView6 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
        
                # -- calc BetaDiversity --
                observeEvent(input$calcBetaDiversity, {
                        
                        rv$infoText_6 <- NULL
                        
                        rv$Tr_PCoA <- NULL
                        rv$adonisTable <- NULL
                        
                        
                        if (is.null(rv$phylum_colors)) {
                                rv$infoText_6 <- "Sorry. You need to calculate phylum colors first."
                                return()
                        }
                        
                        if (input$tcaType6 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType6 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType6 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType6 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_6 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_6 <- "Sorry. The chosen object does not exist yet, did you calculate the size factors?"
                                return()
                        }
                        
                        
                        if (input$filtered6 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_6 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }

                        measure <- input$beta_measure
                        
                        phylum_colors <- rv$phylum_colors
                        
                        # see unlist(phyloseq::distanceMethodList)
                        
                        if (! measure %in% c("jsd", "manhattan", "euclidean", "bray", "jaccard", "canberra")){
                                rv$infoText_6 <- "Sorry. Allowed beta diversity distance measures are: jsd, manhattan, euclidean, bray, caberra, and jaccard today."
                                return()
                                
                        }
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_6 <- error_message
                                return()
                        }
                        # ------
                        

                        # for relative abundance I want it to be relative abundances, for Manhatten it has to be relative abundances otherwise there is an error:
                        
                        if (input$tcaType6 == "library size (RA)" || measure == "manhattan"){
                                ps <- phyloseq::transform_sample_counts(ps, function(x){x/sum(x)})
                        }
                        
                        dist_list <- calc_beta_div_distances(ps, dist_methods = measure, group_var = group_var, compare = group_var_levels)
                        
                        group_factor <- sample_data(ps)[[group_var]]
                        # make sure group_factor fits to dist_list, i.e. only keep samples covered by compare = group_var_levels!
                        group_factor <- factor(group_factor[group_factor %in% group_var_levels], levels = group_var_levels, ordered = T)
                        
                        adonis_list <- lapply(dist_list, function(dist_obj){
                                loop_vegan_adonis(dist_obj = dist_obj, group_fac = group_factor)
                        })
                        
                        
                        if (input$pcoaCorrection == "yes"){
                                coordCor <- TRUE
                        } else {
                                coordCor <- FALSE
                        }
                        
                        pcoas <- calc_ordination_from_distances(ps, group_var = group_var, dist_list = dist_list, color_levels = color_levels, ordination_type = "PCoA", shape = NULL, coord_cor = coordCor, phylum_colors = phylum_colors) 
                        
                        rv$Tr_PCoA <- pcoas[["ordination_Tr_samples"]]
                        
                        rv$adonis <- adonis_list[[1]]
                        
                        rv$infoText_6 <- "Beta diversity analysis has been performed."
                        
                })
                # ----
        
                
                
                # -- render the adonis table --
                output$adonis <- renderTable({
                        if(!is.null(rv$adonis)){
                                adonisShow <- rv$adonis
                                adonisShow
                        } else {
                                NULL
                        }
                }, digits = 4, sanitize.text.function = function(x) x, caption = "adonis test result", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
        
                # -- output plot PCoA --
                output$plotPCoA <- renderPlot({
                        if (is.null(rv$Tr_PCoA)) {
                                NULL
                        } else {
                                rv$Tr_PCoA
                        }
                },
                height = 500)
                # ----
        
        }
        
        ################################ Tab 7 ############################################
        {
                # - Tab7: fisher Test -
                # -- output infoText --
                output$infoText_7 <- renderText({
                        rv$infoText_7
                })
                # ----
                
                
                
                # -- render the overviewViews table --
                output$overviewView7 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
        
                # -- calc fisher test --
                observeEvent(input$fisher, {
                        
                        rv$infoText_7 <- NULL
                        
                        rv$Tr_fisher <- NULL
                        rv$fisherTable <- NULL
                        
                        if (is.null(rv$phylum_colors)) {
                                rv$infoText_7 <- "Sorry. You need to calculate phylum colors first."
                                return()
                        }
                        
                        if (input$tcaType7 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType7 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType7 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType7 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_7 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_7 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                                return()
                        }
                        
                        
                        if (input$filtered7 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_7 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }
                        
                        
                        phylum_colors <- rv$phylum_colors
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_7 <- error_message
                                return()
                        }
                        # ------
                        
                        zero_color <- input$zero_color7
                        
                        if (!areColors(zero_color)){
                                rv$infoText_7 <- "Sorry. The color for the zero counts is not an R color, please change."
                                return()
                                
                        }
                        
                        
                        physeq_to_test <- ps
                        
                        if (input$tcaType7 == "library size (RA)") {
                                physeq_to_test <- phyloseq::transform_sample_counts(physeq_to_test, function(x){x/sum(x)})
                        }
                        
                        
                        diff_ab_df <- test_diffs_in_prevalence_single(physeq = physeq_to_test, group_var = group_var, compare = group_var_levels, p.adj.method = "fdr", minCount = 0L, symnum.args = rv$symnum.args)
                        
                        hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = "fdr")
                        
                        taxa_hit_df <- hit_list[["hit_table"]]
                        
                        significance_colors <- brewer.pal(4, "Reds")
                        significance_colors <- c(rev(significance_colors), "gray", "violet")
                        names(significance_colors) = c("****", "***", "**", "*", "ns", "?")
                        taxa_colors <- list("signi_adj" = significance_colors, "Phylum" = phylum_colors)
                        sample_colors <- list(color_levels)
                        names(sample_colors) <- group_var
                        taxa_annotation <- taxa_hit_df$Annotation
                        
                        
                        max_abundance_for_color <- input$max_abundance_for_colorF
                        
                        if (max_abundance_for_color == ""){
                                max_abundance_for_color <- max(otu_table(physeq_to_test))
                        } else {
                                max_abundance_for_color <- as.numeric(max_abundance_for_color)
                                
                        }
                        
                        
                        if (is.na(max_abundance_for_color) || max_abundance_for_color < min(as(otu_table(physeq_to_test), "matrix")[otu_table(physeq_to_test) > 0]) || max_abundance_for_color > max(otu_table(physeq_to_test), na.rm = TRUE)) {
                                rv$infoText_7 <- "Sorry. max_abundance_for_color either not numeric, too small or too high:). Change it."
                                return()
                                
                        }
                        
                        if (input$log7){
                                logger <- TRUE 
                        } else {
                                logger <- FALSE
                        }
                        
                        max_shown <- as.numeric(input$maxShown7)
                        
                        if (is.na(max_shown) || max_shown < 1) {
                                max_shown <- 40
                        }
                        
                        
                        rv$Tr_fisher <- plot_heatmap_physeq(physeq_to_test, sample_colors = sample_colors, taxa_info_df = head(taxa_hit_df, max_shown), taxa_colors = taxa_colors, 
                                                            taxa_annotation = head(taxa_annotation, max_shown), max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                                                            zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = logger, drop_color_levels = TRUE,
                                                            border_color = NA, 
                                                            cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                                                            annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 22, 
                                                            fontsize_row = 12, fontsize_col = 12, fontsize_number = 12)
                        
                        rv$fisherTable <- taxa_hit_df
                        
                        rv$infoText_7 <- "Fisher test done"
                        
                })
                # ----
        
        
        
                # -- render the fisher result table --
                output$tableFisher <- renderTable({
                        if(!is.null(rv$fisherTable)){
                                fisherShow <- rv$fisherTable
                                fisherShow
                        } else {
                                NULL
                        }
                }, digits = 4, sanitize.text.function = function(x) x, caption = "Results of fisher test result", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
        
                # -- output plot fisher --
                output$plotFisher <- renderPlot({
                        if (is.null(rv$Tr_fisher)) {
                                NULL
                        } else {
                                grid::grid.newpage()
                                grid::grid.draw(rv$Tr_fisher)
                        }
                },
                height = 750)
                # ----
                # --
        }
        
        
        ################################ Tab 8 ############################################
        {
                # - Tab8: DESeq2 Test -
                # -- output infoText --
                output$infoText_8 <- renderText({
                        rv$infoText_8
                })
                # ----
        
                
                # -- render the overviewViews table --
                output$overviewView8 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
        
                # -- calc DESeq2 --
                observeEvent(input$DESeq2, {
                        
                        rv$infoText_8 <- NULL
                        
                        rv$Tr_DESeq2 <- NULL
                        rv$DESeq2Table <- NULL
                        
                        
                        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                                rv$infoText_8 <- "Sorry. You need a loaded raw count phyloseq object, and you need phylum_colors. Make sure these steps are done."
                                return()
                        }
                        
                        ps <- rv$ps
                        
                        if (input$tcaType8 == "DESeq_gm_exclZero") {
                                SFs <- rv$SFs
                        } else if (input$tcaType8 == "DESeq_poscounts") {
                                SFs <- rv$SFs_poscounts
                        } else if (input$tcaType8 == "library size (RA)"){
                                SFs <- rv$SFs_RA
                        }
                        
                        if (is.null(SFs)){
                                rv$infoText_8 <- "Sorry. You need to calculate size factors first!"
                                return()
                                
                        }
                        
                        
                        if (input$filtered8 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_8 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }
                        
                        
                        
                        phylum_colors <- rv$phylum_colors
                        
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_8 <- error_message
                                return()
                        }
                        # ------
                        
                        physeq_to_test <- ps

                        res_list <- test_differential_abundance_DESeq2single(physeq = physeq_to_test, group_var = group_var, compare = group_var_levels, SFs = SFs, p.adjust.method = "fdr", symnum.args = rv$symnum.args, cooksCutoff = TRUE)
                        
                        diff_ab_df <- res_list[[1]]
                        physeq_to_test <- res_list[[2]]
                        
                        hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = "fdr")
                        
                        taxa_hit_df <- hit_list[["hit_table"]]
                        
                        significance_colors <- brewer.pal(4, "Reds")
                        significance_colors <- c(rev(significance_colors), "gray", "violet")
                        names(significance_colors) = c("****", "***", "**", "*", "ns", "?")
                        taxa_colors <- list("signi_adj" = significance_colors, "Phylum" = phylum_colors)
                        sample_colors <- list(color_levels)
                        names(sample_colors) <- group_var
                        taxa_annotation <- taxa_hit_df$Annotation
                        
                        
                        max_abundance_for_color <- input$max_abundance_for_colorD
                        
                        if (max_abundance_for_color == ""){
                                max_abundance_for_color <- max(otu_table(physeq_to_test))
                        } else {
                                max_abundance_for_color <- as.numeric(max_abundance_for_color)
                                
                        }
                        
                        if (is.na(max_abundance_for_color) || max_abundance_for_color < min(as(otu_table(physeq_to_test), "matrix")[otu_table(physeq_to_test) > 0]) || max_abundance_for_color > max(otu_table(physeq_to_test), na.rm = TRUE)) {
                                rv$infoText_8 <- "Sorry. max_abundance_for_color either not numeric, too small or too high:). Change it."
                                return()
                                
                        }
                        
                        zero_color <- input$zero_color8
                        
                        if (!areColors(zero_color)){
                                rv$infoText_8 <- "Sorry. The color for the zero counts is not an R color, please change."
                                return()
                                
                        }
                        
                        
                        if (input$log8){
                                logger <- TRUE 
                        } else {
                                logger <- FALSE
                        }
                        
                        max_shown <- as.numeric(input$maxShown8)
                        
                        if (is.na(max_shown) || max_shown < 1) {
                                max_shown <- 40
                        }
                        
                        rv$Tr_DESeq2 <- plot_heatmap_physeq(physeq_to_test, sample_colors = sample_colors, taxa_info_df = head(taxa_hit_df, max_shown), taxa_colors = taxa_colors, 
                                                            taxa_annotation = head(taxa_annotation, max_shown), max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                                                            zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = logger, drop_color_levels = TRUE,
                                                            border_color = NA, 
                                                            cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                                                            annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 22, 
                                                            fontsize_row = 12, fontsize_col = 12, fontsize_number = 12)
                        
                        rv$DESeq2Table <- taxa_hit_df
                        
                        
                        rv$infoText_8 <- "DESeq2 test done."
                        
                })
                # ----
        
        
        
                # -- render the DESeq2 result table --
                output$tableDESeq2 <- renderTable({
                        if(!is.null(rv$DESeq2Table)){
                                DESeq2Show <- rv$DESeq2Table
                                DESeq2Show
                        } else {
                                NULL
                        }
                }, digits = 4, sanitize.text.function = function(x) x, caption = "Results of DESeq2 analysis", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
        
                # -- output plot DESeq2 --
                output$plotDESeq2 <- renderPlot({
                        if (is.null(rv$Tr_DESeq2)) {
                                NULL
                        } else {
                                grid::grid.newpage()
                                grid::grid.draw(rv$Tr_DESeq2)
                        }
                },
                height = 750)
                # ----
                # --    
        }
        
        
        ################################ Tab 9 ############################################
        {
                # - Tab9: wilcoxon Test -
                # -- output infoText --
                output$infoText_9 <- renderText({
                        rv$infoText_9
                })
                # ----
        
                
                # -- render the overviewViews table --
                output$overviewView9 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
        
                # -- calc wilcoxon --
                observeEvent(input$wilcoxon, {
                        
                        rv$infoText_9 <- NULL
                        
                        rv$Tr_wilcoxon <- NULL
                        rv$wilcoxonTable <- NULL
                        
                        
                        if (is.null(rv$phylum_colors)) {
                                rv$infoText_9 <- "Sorry. You need to calculate phylum colors first."
                                return()
                        }
                        
                        if (input$tcaType9 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType9 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType9 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType9 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_9 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_9 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                                return()
                        }
                        
                        
                        if (input$filtered9 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_9 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }
                        
                        
                        zeros_include <- input$wilcoxonZeros
                        
                        if (zeros_include == "include zero counts") {
                                excludeZeros <- FALSE
                        } else {
                                excludeZeros <- TRUE
                        }
                        
                        
                        phylum_colors <- rv$phylum_colors
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_9 <- error_message
                                return()
                        }
                        # ------
                        
                        physeq_to_test <- ps
                        
                        # physeq_to_test <- phyloseq::transform_sample_counts(physeq_to_test, function(x){x/sum(x)})
                        
                        diff_ab_df <- test_differential_abundance_Wilcoxonsingle(physeq = physeq_to_test, group_var, compare = group_var_levels, 
                                                                                 excludeZeros = excludeZeros, p.adjust.method = "fdr", symnum.args = rv$symnum.args)
                        
                        hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = "fdr")
                        
                        taxa_hit_df <- hit_list[["hit_table"]]
                        
                        significance_colors <- brewer.pal(4, "Reds")
                        significance_colors <- c(rev(significance_colors), "gray", "violet")
                        names(significance_colors) = c("****", "***", "**", "*", "ns", "?")
                        taxa_colors <- list("signi_adj" = significance_colors, "Phylum" = phylum_colors)
                        sample_colors <- list(color_levels)
                        names(sample_colors) <- group_var
                        taxa_annotation <- taxa_hit_df$Annotation
                        
                        max_abundance_for_color <- input$max_abundance_for_colorW
                        
                        if (max_abundance_for_color == ""){
                                max_abundance_for_color <- max(otu_table(physeq_to_test))
                        } else {
                                max_abundance_for_color <- as.numeric(max_abundance_for_color)
                                
                        }
                        
                        if (is.na(max_abundance_for_color) || max_abundance_for_color < min(as(otu_table(physeq_to_test), "matrix")[otu_table(physeq_to_test) > 0]) || max_abundance_for_color > max(otu_table(physeq_to_test), na.rm = TRUE)) {
                                rv$infoText_9 <- "Sorry. max_abundance_for_color either not numeric, too small or too high:). Change it."
                                return()
                                
                        }
                        
                        zero_color <- input$zero_color9
                        
                        if (!areColors(zero_color)){
                                rv$infoText_9 <- "Sorry. The color for the zero counts is not an R color, please change."
                                return()
                                
                        }
                        
                        if (input$log9){
                                logger <- TRUE 
                        } else {
                                logger <- FALSE
                        }
                        
                        max_shown <- as.numeric(input$maxShown9)
                        
                        if (is.na(max_shown) || max_shown < 1) {
                                max_shown <- 40
                        }
                        
                        
                        rv$Tr_wilcoxon <- plot_heatmap_physeq(physeq_to_test, sample_colors = sample_colors, taxa_info_df = head(taxa_hit_df, max_shown), taxa_colors = taxa_colors, 
                                                              taxa_annotation = head(taxa_annotation, max_shown), max_abundance_for_color = max_abundance_for_color, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                                                              zero_color = zero_color, color_function = viridis, color_steps_bw_markers = 10, log_transform = logger, drop_color_levels = TRUE,
                                                              border_color = NA, 
                                                              cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                                                              annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 22, 
                                                              fontsize_row = 12, fontsize_col = 12, fontsize_number = 12)
                        
                        rv$wilcoxonTable <- taxa_hit_df
                        
                        
                        rv$infoText_9 <- "Wilcoxon test done."
                        
                })
                # ----
        
        
        
                # -- render the wilcoxon result table --
                output$tableWilcoxon <- renderTable({
                        if(!is.null(rv$wilcoxonTable)){
                                wilcoxonShow <- rv$wilcoxonTable
                                wilcoxonShow
                        } else {
                                NULL
                        }
                }, digits = 4, sanitize.text.function = function(x) x, caption = "Results of wilcoxon analysis", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
        
                # -- output plot wilcoxon --
                output$plotWilcoxon <- renderPlot({
                        if (is.null(rv$Tr_wilcoxon)) {
                                NULL
                        } else {
                                grid::grid.newpage()
                                grid::grid.draw(rv$Tr_wilcoxon)
                        }
                },
                height = 750)
                # ----
        # --
        }
        
        
        ################################ Tab 10 ############################################
        {
                # - Tab10: Phylum analysis -
                # -- output infoText --
                output$infoText_10 <- renderText({
                        rv$infoText_10
                })
                # ----
                
                # -- render the overviewViews table --
                output$overviewView10 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
                # -- calc Firmicutes to other phyla ratios --
                observeEvent(input$firmicutes, {
                        
                        rv$infoText_10 <- NULL
                        
                        rv$Tr_phylum <- NULL
                        
                        rv$firmicutesTable <- NULL
                        
                        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                                rv$infoText_10 <- "Sorry. A phyloseq raw count object must be present and phylum_colors have to be determined. Do these steps first."
                                return()
                        }
                        
                        if (input$tcaType10 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else {
                                rv$infoText_10 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        
                        if (input$filtered10 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_10 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }
                        
                        
                        ps_phylum <- phyloseq::tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
                        
                        taxonomic_level <- "Phylum"

                        phylum_colors <- rv$phylum_colors
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_10 <- error_message
                                return()
                        }
                        # ------
                        
                        
                        
                        df <- as.data.frame(as(tax_table(ps_phylum), "matrix"))
                        taxa_annotation <- get_taxon_names(df[, 2:7])
                        # taxa_annotation <- strsplit(taxa_annotation, split = "/")
                        # taxa_annotation <- sapply(taxa_annotation, `[`, 1)
                        taxa_annotation <- make.unique(taxa_annotation)
                        
                        taxa_order <- c(names(phylum_colors), taxa_annotation[!taxa_annotation %in% names(phylum_colors)])
                        
                        
                        numerator_phylum <- input$numerator_phylum
                        
                        if (! numerator_phylum %in% df$Phylum){
                                rv$infoText_10 <- "Sorry, numerator_phylum was not a phylum in the phyloseq object, please change the input."
                                return()
                        }
                        
                        
                        denominator_phylum <- input$denominator_phylum
                        
                        if (denominator_phylum == ""){
                                denominator_phylum <- NULL
                        } else {
                                denominator_phylum <- strsplit(x = denominator_phylum, split = ",")
                                denominator_phylum <- unlist(denominator_phylum)
                                denominator_phylum <- gsub(" ", "", denominator_phylum)
                                
                                if (!all(denominator_phylum %in% df$Phylum)){
                                        rv$infoText_10 <- "Sorry, not all given denominator_phyla were a phylum in the phyloseq object, please change the input."
                                        return()
                                }
                        }
                        
                        
                        # NB: if you would like the order based on pValues/significance choose tax_order = NULL
                        FirmicutesRatioPlots <- plot_taxa_ratios_AllLevels(physeq = ps_phylum, group_var = group_var, color_levels = color_levels, tax_names = taxa_annotation, taxa_nom = numerator_phylum, taxa_den = denominator_phylum, tax_order = taxa_order, 
                                                                           test = "wilcox.test", p_adjust_method = "fdr", symnum.args = rv$symnum.args)
                        
                        
                        rv$Tr_phylum <- FirmicutesRatioPlots[["Tr3"]]
                        
                        rv$firmicutesTable <- FirmicutesRatioPlots[["pValsLog"]]
                        
                        
                        rv$infoText_10 <- "Ratios of Firmicutes to all other phyla have been calculated."
                        
                })
                # ----
        
        
                # -- calc phylum to phylum tile plots --
                observeEvent(input$tile, {
                        
                        rv$infoText_10 <- NULL
                        
                        rv$Tr_phylum <- NULL
                        
                        rv$firmicutesTable <- NULL
                        
                        if (is.null(rv$phylum_colors) || is.null(rv$ps)) {
                                rv$infoText_10 <- "Sorry. A phyloseq raw count object must be present and phylum_colors have to be determined. Do these steps first."
                                return()
                        }
                        
                        if (input$tcaType10 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else {
                                rv$infoText_10 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        
                        if (input$filtered10 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_10 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }
                        
                        
                        ps_phylum <- phyloseq::tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
                        
                        
                        taxonomic_level <- "Phylum"
                        
                        
                        phylum_colors <- rv$phylum_colors
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_10 <- error_message
                                return()
                        }
                        # ------
                        
                        
                        
                        raw_TbTmatrixes <- calculate_raw_TbTmatrixes(physeq = ps_phylum)
                        
                        
                        df <- as.data.frame(as(tax_table(ps_phylum), "matrix"))
                        taxa_annotation <- get_taxon_names(df[, 2:7])
                        # taxa_annotation <- strsplit(taxa_annotation, split = "/")
                        # taxa_annotation <- sapply(taxa_annotation, `[`, 1)
                        taxa_annotation <- make.unique(taxa_annotation)
                        
                        taxa_order <- c(names(phylum_colors), taxa_annotation[!taxa_annotation %in% names(phylum_colors)])
                        
                        TbT_tile <- create_raw_TbT_TilePlot(TbTmatrixes = raw_TbTmatrixes, physeq = ps_phylum, group_var = group_var, color_levels = color_levels, signi_level = 0.05, tax_names = taxa_annotation, tax_order = taxa_order, test = "wilcoxon", p_adjust_method = "none")
                        
                        rv$Tr_phylum <- TbT_tile
                        
                        
                        rv$infoText_10 <- "Tile plot of phylum/phylum ratios has been generated. If a tile is colored with the color of a sample group, this group has higher (row phylum)/(column phylum) ratios."
                        
                })
                # ----
        
        
        
                # -- render the Firmicutes result table --
                output$tableFirmicutes <- renderTable({
                        if(!is.null(rv$firmicutesTable)){
                                firmicutesShow <- rv$firmicutesTable
                                firmicutesShow
                        } else {
                                NULL
                        }
                }, digits = 5, sanitize.text.function = function(x) x, caption = "Results of Firmicutes ratios analysis", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
        
        
        
                # -- output plot phylum analysis --
                output$plotPhylum <- renderPlot({
                        if (is.null(rv$Tr_phylum)) {
                                rv$Tr_phylum
                        } else {
                                rv$Tr_phylum
                        }
                },
                height = 750)
                # ----
        # --
        }
        ################################ Tab 13 ############################################
        {
                # - Tab13: Taxa finder -
                # -- output infoText --
                output$infoText_13 <- renderText({
                        rv$infoText_13
                })
                # ----
                
                # -- render the overviewViews table --
                output$overviewView13 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                # -- search Taxa --
                observeEvent(input$taxaSearch, {
                        
                        rv$infoText_13 <- NULL
                        
                        rv$taxaFindTable <- NULL
                        
                        rv$Tr_taxaAbundance <- NULL
                        
                        rv$tablepValsTaxa <- NULL
                        
                        
                        if (input$tcaType13 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType13 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType13 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType13 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_13 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_13 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                                return()
                        }
                        
                        
                        if (input$filtered13 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_13 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }
                        
                        tax_level <- input$taxLevelSearch
                        
                        TT <- as.data.frame(unclass(tax_table(ps)))
                        
                        searchWord <- input$searchWord
                        
                        indexes <- grep(pattern = searchWord, TT[[tax_level]], ignore.case = TRUE)
                        
                        if (length(indexes) == 0){
                                rv$infoText_13 <- "No taxon matched your search."
                                return()
                        }
                        
                        TT_show <- TT[indexes,]
                        
                        TT_show <- cbind(data.frame(Index = indexes, Annotation = get_taxon_names(TT_show)), TT_show)
                        
                        rv$taxaFindTable <- TT_show
                        
                        rv$infoText_13 <- "Taxa have been found"

                        # ------
                        
                })
                # ----
                
                
                
                # -- search Taxa --
                observeEvent(input$plotTaxa, {
                        
                        rv$infoText_13 <- NULL
                        
                        rv$Tr_taxaAbundance <- NULL
                        
                        rv$tablepValsTaxa <- NULL
                        
                        if (is.null(rv$phylum_colors)) {
                                rv$infoText_13 <- "Sorry. You need to calculate phylum colors first."
                                return()
                        }
                        
                        
                        if (input$tcaType13 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType13 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType13 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType13 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_13 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_13 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                                return()
                        }
                        
                        
                        if (input$filtered13 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_13 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }
                        
                
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_13 <- error_message
                                return()
                        }
                        # ------
                        
                        
                        userTaxa <- input$userTaxa
                        
                        userTaxa <- strsplit(x = userTaxa, split = ",")
                        
                        userTaxa <- as.numeric(unlist(userTaxa))
                        
                        userTaxa <- userTaxa[!is.na(userTaxa)]
                        
                        if (length(userTaxa) == 0){
                                rv$infoText_13 <- "No correct indexes were specified. Please give a comma-separated list of the indexes of the taxa you want to plot."
                                return()
                        }
                        
                        userTaxa <- sort(userTaxa)
                        
                        if (!all(userTaxa %in% 1:ntaxa(ps))){
                                rv$infoText_13 <- "Some indexes were outside of the range of the chosen phyloseq object so please correct your input and try again."
                                return()
                        }
                        
                        keepTaxa <- taxa_names(ps)[userTaxa]
                        
                        ps.pruned <- phyloseq::prune_taxa(taxa = keepTaxa, ps)

                        mdf <- psmelt(ps.pruned)
                        
                        TT <- as.data.frame(unclass(tax_table(ps.pruned)))
                        TT$Annotation <- get_taxon_names_plusTL(TT)
                        LookUp <- data.frame(Index = userTaxa, OTU = keepTaxa, Annotation = TT$Annotation[match(keepTaxa, rownames(TT))])
                        
                        LookUp <- dplyr::mutate(LookUp, Labeller = paste(Index, OTU, Annotation, sep = "_"))
                        
                        mdf$Index <- LookUp$Index[match(mdf$OTU, LookUp$OTU)]
                        mdf$Annotation <- LookUp$Annotation[match(mdf$OTU, LookUp$OTU)]
                        mdf$Labeller <- paste(mdf$Index, mdf$OTU, mdf$Annotation, sep = "_")
                        
                        mdf$Labeller <- factor(mdf$Labeller, levels = LookUp$Labeller, ordered = TRUE) 
                        
                        mdf <- mdf[mdf[[group_var]] %in% group_var_levels,]
                        
                        if (input$pool){
                                
                                mdf$Labeller <- paste0("Indexes_", paste(userTaxa, collapse = "_"))
                                
                                mdf <- group_by_(mdf, "Sample", "Phylum", group_var, "Labeller")
                                
                                mdf <- summarise(mdf, Abundance = sum(Abundance))
                        }
                        
                        
                        
                        mdf[[group_var]] <- factor(mdf[[group_var]], levels = names(color_levels), ordered = TRUE)
                        
                        group_fac <- factor(mdf[[group_var]])
                        fac_levels <- levels(group_fac)
                        
                        comparisonList <- get_unique_facLevel_combinations(fac_levels)
                        
                        phylum_colors <- rv$phylum_colors
                        
                        Tr <- ggplot(mdf, aes_string(x = group_var, y = "Abundance", fill = "Phylum", col = group_var))
                        
                        Tr <- Tr +
                                geom_boxplot(outlier.colour = NA) +
                                geom_jitter(heigth = 0) +
                                facet_wrap(~ Labeller, scales = "free_y") +
                                scale_color_manual(values = color_levels) +
                                scale_fill_manual(values = phylum_colors) +
                                theme_bw() +
                                xlab("")
                        
                        if (input$log){
                                Tr <- Tr + scale_y_log10()
                                mdf$Abundance <- log(mdf$Abundance)
                                mdf <- mdf[is.finite(mdf$Abundance),]
                        }
                        
                        
                        
                        if (all(table(mdf[[group_var]]) > 2)) {
                                
                                Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = "t.test", hide.ns = FALSE, symnum.args = rv$symnum.args)
                                
                                formulaa <- as.formula(paste("Abundance ~", group_var, sep = " "))
                                
                                pVals <- compare_means(formula = formulaa, data = mdf, group.by = "Labeller", method = "t.test", p.adjust.method = "fdr", symnum.args = rv$symnum.args)
                                
                        } else {
                                pVals <- NULL
                        }
                        
                        rv$Tr_taxaAbundance <- Tr
                        
                        rv$tablepValsTaxa <- pVals
                        
                        rv$infoText_13 <- "selected Taxa have been plotted."
                        
                        # ------
                        
                })
                # ----
                
                
                
                
                # -- render the foundTaxa table --
                output$foundTaxa <- renderTable({
                        if(!is.null(rv$taxaFindTable)){
                                taxaShow <- rv$taxaFindTable
                                taxaShow
                        } else {
                                NULL
                        }
                }, rownames = TRUE, digits = 0, sanitize.text.function = function(x) x, caption = "Taxa found in your search", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                # -- render the pVals of the taxa comparison --
                output$pValsTaxa <- renderTable({
                        if(!is.null(rv$tablepValsTaxa)){
                                rv$tablepValsTaxa
                        } else {
                                NULL
                        }
                }, rownames = TRUE, digits = 0, sanitize.text.function = function(x) x, caption = "p Value table", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                
                # -- output plot Taxa Abundances --
                output$plotTaxaAbundances <- renderPlot({
                        if (is.null(rv$Tr_taxaAbundance)) {
                                rv$Tr_taxaAbundance
                        } else {
                                rv$Tr_taxaAbundance
                        }
                }, 
                height = 700)
                # ----
                
        }
        
        ################################ Tab 14 ############################################
        {
                # - Tab14: Taxa ratio -
                # -- output infoText --
                output$infoText_14 <- renderText({
                        rv$infoText_14
                })
                # ----
                
                
                # -- render the overviewViews table --
                output$overviewView14 <- renderTable({
                        rv$overviewView
                }, sanitize.text.function = function(x) x, caption = "Overview of phyloseq objects and filtered taxa", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                
                # -- plot Ratio Taxa --
                observeEvent(input$plotRatioTaxa, {
                        
                        rv$infoText_14 <- NULL
                        
                        rv$Tr_taxaRatios <- NULL
                        
                        rv$tablepValsTaxaRatios <- NULL
                        
                        if (input$tcaType14 == "raw counts"){
                                ps <- rv$ps
                                # sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]], levels = c(group_var_levels, setdiff(unique(sample_data(ps)[[group_var]]), group_var_levels)), ordered = TRUE)
                                # rv$ps <- ps
                        } else if (input$tcaType14 == "DESeq_gm_exclZero"){
                                ps <- rv$ps_tca
                        } else if (input$tcaType14 == "DESeq_poscounts") {
                                ps <- rv$ps_tca_poscounts
                        } else if (input$tcaType14 == "library size (RA)") {
                                ps <- rv$ps_RA
                        } else {
                                rv$infoText_14 <- "Weird how could you choose non of the given options?"
                                return()
                        }
                        
                        if (is.null(ps)) {
                                rv$infoText_14 <- "Sorry. The chosen object does not exist yet. Did you load one?"
                                return()
                        }
                        
                        
                        if (input$filtered14 == "filtered taxa") {
                                
                                if (is.null(rv$filteredKeepTaxa)){
                                        rv$infoText_14 <- "Sorry. You have to do the filtering first."
                                        return()
                                }
                                
                                keepTaxa <- rv$filteredKeepTaxa
                                ps <- prune_taxa(keepTaxa, ps)
                        }
                        
                        
                        # --- the input parameter check ---
                        group_var <- rv$group_var
                        
                        group_var_levels <- rv$group_var_levels 
                        
                        color_levels <- rv$color_levels
                        
                        error_message <- check_user_parameters(group_var = group_var, group_var_levels = group_var_levels, 
                                                               color_levels = color_levels, ps = ps)
                        
                        if (!is.null(error_message)){
                                rv$infoText_14 <- error_message
                                return()
                        }
                        # ------
                        
                        
                        numerator <- as.numeric(input$numerator)
                        
                        if (length(numerator) == 0 || is.na(numerator) || !(numerator %in% 1:ntaxa(ps))){
                                rv$infoText_14 <- "Couldn't find numerator. Change input."
                                return()
                        }
                        
                        
                        denominator <- as.numeric(input$denominator)
                        
                        if (length(denominator) == 0 || is.na(denominator) || !(denominator %in% 1:ntaxa(ps))){
                                rv$infoText_14 <- "Couldn't find denominator. Change input."
                                return()
                        }
                        
                        
                        keepTaxa <- taxa_names(ps)[c(numerator, denominator)]
                        
                        ps.pruned <- phyloseq::prune_taxa(taxa = keepTaxa, ps)
                        
                        mdf <- psmelt(ps.pruned)
                        
                        TT <- as.data.frame(unclass(tax_table(ps.pruned)))
                        TT$Annotation <- get_taxon_names_plusTL(TT)
                        LookUp <- data.frame(Index = c(numerator, denominator), OTU = keepTaxa, Annotation = TT$Annotation[match(keepTaxa, rownames(TT))])
                        
                        mdf$Index <- LookUp$Index[match(mdf$OTU, LookUp$OTU)]
                        
                        mdf <- mdf[mdf[[group_var]] %in% group_var_levels,]
                        
                        mdf_ratio <- group_by_(mdf, "Sample", group_var)
                        
                        mdf_ratio <- summarise(mdf_ratio, Ratio = Abundance[Index == numerator]/Abundance[Index == denominator])
                        
                        mdf_ratio$Labeller <- paste0(numerator, "/", denominator, "_", paste(LookUp$Annotation, collapse = "/"))
                        
                        mdf_ratio <- mdf_ratio[is.finite(mdf_ratio$Ratio),]
                        
                        
                        mdf_ratio[[group_var]] <- factor(mdf_ratio[[group_var]], levels = names(color_levels), ordered = TRUE)
                        
                        group_fac <- factor(mdf_ratio[[group_var]])
                        fac_levels <- levels(group_fac)
                        
                        comparisonList <- get_unique_facLevel_combinations(fac_levels)
                        
                        Tr <- ggplot(mdf_ratio, aes_string(x = group_var, y = "Ratio", col = group_var))
                        
                        Tr <- Tr +
                                geom_boxplot(outlier.colour = NA) +
                                geom_jitter(height = 0) +
                                facet_wrap(~ Labeller, scales = "free_y") +
                                scale_color_manual(values = color_levels) +
                                theme_bw() +
                                xlab("")
                        
                        if (input$logRatio){
                                Tr <- Tr + scale_y_log10()
                                mdf_ratio$Ratio <- log(mdf_ratio$Ratio)
                                mdf_ratio <- mdf_ratio[is.finite(mdf_ratio$Ratio),]
                        }
                        
                        
                        
                        
                        if (all(table(mdf_ratio[[group_var]]) > 2)) {
  
                                Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = "t.test", hide.ns = FALSE, symnum.args = rv$symnum.args)
                                
                                formulaa <- as.formula(paste("Ratio ~", group_var, sep = " "))
                                
                                pVals <- compare_means(formula = formulaa, data = mdf_ratio, group.by = "Labeller", method = "t.test", p.adjust.method = "fdr", symnum.args = rv$symnum.args)
                        } else {
                                pVals <- NULL
                        }
                        
                        rv$Tr_taxaRatios <- Tr
                        
                        rv$tablepValsTaxaRatios <- pVals
                        
                        
                        rv$infoText_14 <- "Selected abundance ratios have been plotted and statistical test has been performed if possible."
                        
                        # ------
                        
                })
                # ----
                
                
                
                # -- render the pVals of the ratio comparison --
                output$pValsRatioTaxa <- renderTable({
                        if(!is.null(rv$tablepValsTaxaRatios)){
                                rv$tablepValsTaxaRatios
                        } else {
                                NULL
                        }
                }, rownames = TRUE, digits = 0, sanitize.text.function = function(x) x, caption = "p Value table", caption.placement = getOption("xtable.caption.placement", "top"))
                # ----
                
                
                
                # -- output plot Taxa Abundances --
                output$plotTaxaRatios <- renderPlot({
                        if (is.null(rv$Tr_taxaRatios)) {
                                rv$Tr_taxaRatios
                        } else {
                                rv$Tr_taxaRatios
                        }
                }, 
                height = 500)
                # ----
                
        }
        
        
}

## ======== Run the app ===================

shinyApp(ui = ui, server = server)
