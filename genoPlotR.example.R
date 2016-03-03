## genoPlotR.example.R

## code to run genoPlotR (comparative gene and genome visualization) and plot CBET3_02473 cluster alignment
##installation:
install.packages("genoPlotR", repos="http://R-Forge.R-project.org")

##loading the package
library(genoPlotR)

setwd("../analyses/antismash2/CBET3_02473")

#load the dnas (annotation files)
CBE1 <- read_dna_seg_from_tab(("CBE2473.annots"), header = TRUE)
CBE2 <- read_dna_seg_from_tab(("CBE2338.annots"), header = TRUE)
PTR1 <- read_dna_seg_from_tab(("PTR2708.annots"), header = TRUE)
ANI1 <- read_dna_seg_from_tab(("ANI6832.annots"), header = TRUE)

#define middle of a feature
mpCBE1 <- middle(CBE1)
mpCBE2 <- middle(CBE2)
mpPTR1 <- middle(PTR1)
mpANI1 <- middle(ANI1)

#add annotations to the middle
annotCBE1 <- annotation(x1=mpCBE1, text=CBE1$name, rot=35)
annotCBE2 <- annotation(x1=mpCBE2, text=CBE2$name, rot=35)
annotPTR1 <- annotation(x1=mpPTR1, text=PTR1$name, rot=35)
annotANI1 <- annotation(x1=mpANI1, text=ANI1$name, rot=35)

#lengths (and orientations) of the clusters
compSegs <- list(c(0,46292), c(0,21849), c(0,43020), c(0,43957))

#import the bl2seq comparisons (direction in bl2seq (defined by -i and -j) is crucial, headers should be clean (no spaces etc.))
CBE1_vs_CBE2 <- try(read_comparison_from_blast("CBE1_vs_CBE2.blast",
                                            filt_length = 100,))
CBE2_vs_PTR1 <- try(read_comparison_from_blast("CBE2_vs_PTR.blast",
                                               filt_length = 100,))
PTR1_vs_ANI1 <- try(read_comparison_from_blast("PTR_vs_ANI.blast",
                                               filt_length = 100,))

#define comparisons, annotations, dnas and names
comparisons <- list(CBE1_vs_CBE2,CBE2_vs_PTR1,PTR1_vs_ANI1)
annotations <- list(annotCBE1,annotCBE2,annotPTR1,annotANI1)
dna_segs <- list(CBE1,CBE2,PTR1,ANI1)
names <- c("CBET3_02473", "CBET3_02338", "PTR2708","ANI6832")
names(dna_segs) <- names

#plot the cluster alignment using genoPlotR
plot_gene_map(dna_segs = dna_segs, 
              comparisons = comparisons, 
              gene_type = "arrows", 
              dna_seg_scale = TRUE, 
              scale = FALSE,
              limit_to_longest_dna_seg = FALSE, 
              main_pos="centre",
              xlims=compSegs,
              annotations=annotations,
              annotation_cex=0.6,
              annotation_height=2.5,
)
