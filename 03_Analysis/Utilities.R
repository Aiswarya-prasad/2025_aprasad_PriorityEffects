working_dir <- "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/"
setwd(working_dir)
# use R in /work/FAC/FBM/DMF/pengel/spirit/aprasad/miniconda3/envs/pacbio_ampli_r_env/bin/R
library(ggplot2)
library(readxl)
library(knitr)
library(RColorBrewer)
library(scales)
library(dplyr)
library(vegan)
library(ggsignif)
library(plotly)
library(htmlwidgets)
library(viridis)
library(hrbrthemes)
library(ggthemes)
library(ggrepel)
library(gridExtra)
library(microbiome)
library(ape)
library(phyloseq)
library(ComplexHeatmap)
library(ggnewscale)
library(VennDiagram)
library(ggforce)
library(circlize)
library(magrittr)
library(tidyverse)
library(corrplot)
# library(ggfx) # needed mamba install -c conda-forge r-systemfonts>=1.1.0 but did not happen
library(igraph)
library(tidygraph)
library(ggdendro)
library(dendextend)
library(ggpubr)
library(gggenes)
library(dada2)
library(ggbeeswarm)
library(lme4)
library(rstatix)
library(Biostrings)


make_theme <- function(theme_name=theme_classic(), max_colors=0, palettefill="Pastel1", palettecolor="Dark2", modify_guide = T,
                        setFill=TRUE, setCol=TRUE,
                        guide_nrow=2, guide_nrow_byrow=TRUE, leg_pos="top", leg_size=12,
                        axis_x_title = 12, axis_y_title = 12,
                        x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                        y_angle=0 ,y_vj=0, y_hj=0, y_size=12){
  n_11 = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  n_12 = c("Paired", "Set3")
  n_8 = c("Accent", "Dark2", "Pastel2", "Set2")
  if (palettefill %in% n_12) {
    n_f = 12
  } else {
    if (palettefill %in% n_11) {
      n_f = 11
    } else {
      if (palettefill %in% n_8) {
        n_f  = 8
      } else {
        n_f = 9
      }
    }
  }
  if (palettecolor %in% n_12) {
    n_c = 12
  } else {
    if (palettecolor %in% n_11) {
      n_c = 11
    } else {
      if (palettecolor %in% n_8) {
        n_c  = 8
      } else {
        n_c = 9
      }
    }
  }
  getFill = colorRampPalette(brewer.pal(n_f, palettefill))
  getColor = colorRampPalette(brewer.pal(n_c, palettecolor))
  theme_params <- theme(axis.text.x = element_text(angle = x_angle,
    vjust = x_vj, hjust=x_hj,
    size = x_size),
    axis.text.y = element_text(angle = y_angle,
      vjust = y_vj, hjust=y_hj,
      size = y_size),
      axis.title.x = element_text(size=axis_x_title),
      axis.title.y = element_text(size=axis_y_title),
      legend.position=leg_pos,
      legend.text = element_text(size=leg_size)
    )
  if (modify_guide == T) {
    guide_params <- guides(fill = guide_legend(
                                    nrow=guide_nrow,
                                    byrow=guide_nrow_byrow
                                  ),
                          col = guide_legend(
                                    nrow=guide_nrow,
                                    byrow=guide_nrow_byrow
                                  )
                    )
  my_theme <- list(
                theme_name,
                theme_params,
                guide_params
              )
  } else {
    my_theme <- list(
                  theme_name,
                  theme_params
                )
  }
  if(setFill) {
    if (n_f < max_colors) {
      my_theme <- list(
                    my_theme,
                    scale_fill_manual(values = getFill(max_colors), na.value="grey")
                  )

    } else {
      my_theme <- list(
                    my_theme,
                    scale_fill_brewer(palette=palettefill, na.value="grey")
                  )
    }
  }
  if(setCol) {
    if (n_c < max_colors) {
      my_theme <- list(
                    my_theme,
                    scale_color_manual(values = getColor(max_colors), na.value="grey")
                  )

    } else {
      my_theme <- list(
                    my_theme,
                    scale_color_brewer(palette=palettecolor, na.value="grey")
                  )
    }
  }
  return(my_theme)
}


sample_name_from_shortname <- function(shortname) {
  sample_name <- gsub("-", "_", shortname)
  return(paste0("DNA", sample_name))
}

strains_all = c("ESL0825","ESL0820",
"ESL0822","ESL0170",
"ESL0824","ESL0198",
"ESL0827","ESL0199",
"ESL0200","ESL0819",
"ESL0197",
"ESL0294","ESL0393","ESL0295",
"ESL1028",
"ESL0353","ESL0263","ESL0185",
"ESL0183","ESL0262","ESL0835","ESL0354",
"ESL0394","ESL0259","ESL0260","ESL0184","ESL0350",
"ESL0186","ESL0261","ESL0351")

# strains handled combines the identical strains so they do not complicate the 
# correctness of formulating this as a lin alg problem (as the strains are indistinguishable)
strains_handled = c("ESL0825","ESL0820",
"ESL0822","ESL0170",
"ESL0824","ESL0198",
"ESL0827","ESL0199",
"ESL0200","ESL0819",
"ESL0197",
"ESL0294","ESL0393","ESL0295",
"ESL1028",
"ESL0353_ESL0185","ESL0263",
"ESL0183_ESL0262","ESL0835","ESL0354",
"ESL0394","ESL0259","ESL0260","ESL0184","ESL0350",
"ESL0186","ESL0261","ESL0351")

remove_extension <- function(x, ext) {
  return(gsub(ext, "", x))
}

strains = c("ESL0825",
"ESL0820",
"ESL0822",
"ESL0170",
"ESL0824",
"ESL0198",
"ESL0827",
"ESL0199",
"ESL0200",
"ESL0819",
"ESL0197",
"ESL0295",
"ESL1028",
"ESL0185",
"ESL0263",
"ESL0183",
"ESL0835",
"ESL0394",
"ESL0184",
"ESL0350",
"ESL0186",
"ESL0261")

strains_uniq = strains
# strains_uniq = c("ESL0825",
# "ESL0820",
# "ESL0822",
# "ESL0170",
# "ESL0824",
# "ESL0198",
# "ESL0827",
# "ESL0199",
# "ESL0200",
# "ESL0819",
# "ESL0197",
# "ESL0295",
# "ESL0394",
# "ESL1028",
# "ESL0353_ESL0185",
# "ESL0185",
# "ESL0263",
# "ESL0183_ESL0262",
# "ESL0183",
# "ESL0835",
# "ESL0184",
# "ESL0350",
# "ESL0186",
# "ESL0261")

strain_without_species_pair <- c("ESL0197",
                                 "ESL1028",
                                 "ESL0351",
                                 "ESL0295",
                                 "unknown")

strainColorsUniq = c(
    # Bifidobacterium
    "ESL0825" = "#72a0c1", # Bifidobacterium apousia
    "ESL0820" = "#003366", # Bifidobacterium apousia
    "ESL0822" = "#a6cee3", # Bifidobacterium asteroides
    "ESL0170" = "#1f78b4", # Bifidobacterium asteroides
    "ESL0824" = "#b2df8a", # Bifidobacterium polysaccharolyticum
    "ESL0198" = "#33a02c", # Bifidobacterium polysaccharolyticum
    "ESL0827" = "#8dd3c7", # Bifidobacterium sp1.
    "ESL0199" = "#009080", # Bifidobacterium sp1.
    "ESL0200" = "#ccebc5", # Bifidobacterium sp2.
    "ESL0819" = "#4d7c4f", # Bifidobacterium sp2.
    "ESL0197" = "#b9b9b9", # Bifidobacterium coryneforme

    # Bombilactobacillus
    "ESL0295" = "#cab2f6", # Bombilactobacillus mellis
    "ESL1028" = "#b9b9b9", # Bombilactobacillus mellifer

    # Lactobacillus
    "ESL0185" = "#ffff99", # Lactobacillus apis
    "ESL0353_ESL0185" = "#ffff99", # Lactobacillus apis
    "ESL0263" = "#b15928", # Lactobacillus apis
    "ESL0183" = "#fcade5", # Lactobacillus helsingborgensis
    "ESL0183_ESL0262" = "#fcade5", # Lactobacillus helsingborgensis
    "ESL0835" = "#e7298a", # Lactobacillus helsingborgensis
    "ESL0184" = "#fb9a99", # Lactobacillus melliventris
    "ESL0350" = "#e31a1c", # Lactobacillus melliventris
    "ESL0394" = "#a50f15", # Lactobacillus melliventris
    
    "ESL0186" = "#fdbf6f", # Lactobacillus kullabergensis
    "ESL0261" = "#ff7f00"  # Lactobacillus kullabergensis
)


strain_name_spec_dict_uniq = c("ESL0825" = "Bifidobacterium apousia-ESL0825",
                        "ESL0820" = "Bifidobacterium apousia-ESL0820",
                        "ESL0822" = "Bifidobacterium asteroides-ESL0822",
                        "ESL0170" = "Bifidobacterium asteroides-ESL0170",
                        "ESL0824" = "Bifidobacterium polysaccharolyticum-ESL0824",
                        "ESL0198" = "Bifidobacterium polysaccharolyticum-ESL0198",
                        "ESL0827" = "Bifidobacterium sp1.-ESL0827",
                        "ESL0199" = "Bifidobacterium sp1.-ESL0199",
                        "ESL0200" = "Bifidobacterium sp2.-ESL0200",
                        "ESL0819" = "Bifidobacterium sp2.-ESL0819",
                        "ESL0197" = "Bifidobacterium coryneforme-ESL0197",
                        "ESL0294" = "Bombilactobacillus mellis-ESL0294",
                        "ESL0295" = "Bombilactobacillus mellis-ESL0295",
                        
                        "ESL0393" = "Bombilactobacillus mellis-ESL0393",

                        "ESL1028" = "Bombilactobacillus mellifer-ESL1028",
                        "ESL0353" = "Lactobacillus apis-ESL0353",
                        "ESL0263" = "Lactobacillus apis-ESL0263",
                        "ESL0185" = "Lactobacillus apis-ESL0185",
                        "ESL0353_ESL0185" = "Lactobacillus apis-ESL0185",
                        "ESL0183" = "Lactobacillus helsingborgensis-ESL0183",
                        "ESL0183_ESL0262" = "Lactobacillus helsingborgensis-ESL0183",
                        "ESL0262" = "Lactobacillus helsingborgensis-ESL0262",
                        "ESL0835" = "Lactobacillus helsingborgensis-ESL0835",
                        "ESL0354" = "Lactobacillus helsingborgensis-ESL0354",
                        "ESL0259" = "Lactobacillus melliventris-ESL0259",
                        "ESL0260" = "Lactobacillus melliventris-ESL0260",
                        "ESL0184" = "Lactobacillus melliventris-ESL0184",
                        "ESL0350" = "Lactobacillus melliventris-ESL0350",
                        
                        "ESL0394" = "Lactobacillus melliventris-ESL0394",
                        
                        "ESL0186" = "Lactobacillus kullabergensis-ESL0186",
                        "ESL0261" = "Lactobacillus kullabergensis-ESL0261",
                        "ESL0351" = "Lactobacillus kullabergensis-ESL0351")

species_name_spec_dict = c("ESL0825" = "Bifidobacterium apousia",
                        "ESL0820" = "Bifidobacterium apousia",
                        "ESL0822" = "Bifidobacterium asteroides",
                        "ESL0170" = "Bifidobacterium asteroides",
                        "ESL0824" = "Bifidobacterium polysaccharolyticum",
                        "ESL0198" = "Bifidobacterium polysaccharolyticum",
                        "ESL0827" = "Bifidobacterium sp1.",
                        "ESL0199" = "Bifidobacterium sp1.",
                        "ESL0200" = "Bifidobacterium sp2.",
                        "ESL0819" = "Bifidobacterium sp2.",
                        "ESL0197" = "Bifidobacterium coryneforme",
                        "ESL0294" = "Bombilactobacillus mellis",
                        "ESL0295" = "Bombilactobacillus mellis",
                        
                        "ESL0393" = "Bombilactobacillus mellis",
                        
                        "ESL1028" = "Bombilactobacillus mellifer",
                        "ESL0353" = "Lactobacillus apis",
                        "ESL0263" = "Lactobacillus apis",
                        "ESL0185" = "Lactobacillus apis",
                        "ESL0353_ESL0185" = "Lactobacillus apis",
                        "ESL0183" = "Lactobacillus helsingborgensis",
                        "ESL0183_ESL0262" = "Lactobacillus helsingborgensis",
                        "ESL0262" = "Lactobacillus helsingborgensis",
                        "ESL0835" = "Lactobacillus helsingborgensis",
                        "ESL0354" = "Lactobacillus helsingborgensis",
                        "ESL0259" = "Lactobacillus melliventris",
                        "ESL0260" = "Lactobacillus melliventris",
                        "ESL0184" = "Lactobacillus melliventris",
                        "ESL0350" = "Lactobacillus melliventris",
                        
                        "ESL0394" = "Lactobacillus melliventris",
                        
                        "ESL0186" = "Lactobacillus kullabergensis",
                        "ESL0261" = "Lactobacillus kullabergensis",
                        "ESL0351" = "Lactobacillus kullabergensis")

genusColors <- c("Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                 "Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                 "Bifidobacterium" = brewer.pal(11, "Spectral")[3]
)

speciesColors_uniq <- c("Bifidobacterium apousia" = "#003366",
                  "Bifidobacterium asteroides" = "#1f78b4",
                  "Bifidobacterium polysaccharolyticum" = "#33a02c",
                  "Bifidobacterium sp1." = "#009080",
                  "Bifidobacterium sp2." = "#4d7c4f",
                  "Bifidobacterium coryneforme" = "#c7e9c0",
                  "Bombilactobacillus mellis" = "#6a3d9a",
                  "Bombilactobacillus mellifer" = "#dadaeb",
                  "Lactobacillus apis" = "#b15928",
                  "Lactobacillus helsingborgensis" = "#e7298a",
                  "Lactobacillus melliventris" = "#e31a1c",
                  "Lactobacillus kullabergensis" = "#ff7f00",
                  "unknown" = "#999999")

# for each species its strains get dark and light of the same shade

strainColors_ext <- c("ESL0825" = "#6a3d9a",
                   "ESL0820" = "#cab2d6",
                   "ESL0822" = "#6a3d9a",
                   "ESL0170" = "#cab2d6",
                   "ESL0824" = "#6a3d9a",
                   "ESL0198" = "#cab2d6",
                   "ESL0827" = "#6a3d9a",
                   "ESL0199" = "#cab2d6",
                   "ESL0200" = "#6a3d9a",
                   "ESL0819" = "#cab2d6",
                   "ESL0197" = "#999999",
                   "ESL0294" = "#33a02c",
                   "ESL0393" = "#b2df8a",
                   "ESL0295" = "#00441b",
                   "ESL1028" = "#999999",
                   "ESL0353_ESL0185" = "#e31a1c",
                   "ESL0185" = "#e31a1c",
                   "ESL0353" = "#e31a1c",
                   "ESL0263" = "#fb9a99",
                   "ESL0835" = "#e31a1c",
                   "ESL0183_ESL0262" = "#fb9a99",
                   "ESL0183" = "#fb9a99",
                   "ESL0262" = "#fb9a99",
                   "ESL0354" = "#ff7f00",
                   "ESL0394" = "#e31a1c",
                   "ESL0259" = "#fb9a99",
                   "ESL0260" = "#ff7f00",
                   "ESL0184" = "#fdbf6f",
                   "ESL0350" = "#a50f15",
                   "ESL0186" = "#e31a1c",
                   "ESL0261" = "#fb9a99",
                   "ESL0351" = "#ff7f00",
                   "unknown" = "#999999")

all_treatments <- c(
    "MD", "extraction_blank", "A1-0",
    "B1-0", "A1-B1", "B1-A1",
    "A2-0", "B2-0", "A2-B2",
    "B2-A2", "A3-0", "B3-0",
    "A3-B3", "B3-A3", "A1",
    "A2", "A3", "B1",
    "B2", "B3", "SW",
    "A4-0", "B4-0", "A4-B4",
    "B4-A4", "A5-0", "B5-0",
    "A5-B5", "B5-A5", "A6-0",
    "B6-0", "A6-B6", "B6-A6",
    "A4", "A5", "A6",
    "B4", "B5", "B6",
    "Am1-B1", "Bm1-A1", "Am4-B1",
    "Bm4-A1", "Am6-B1", "Bm6-A1",
    "dropout_community", "Am7-B1", "Bm7-A1"
)

all_exp_treatments <- c(
    "MD",
    "A1-0", "B1-0", "A1-B1", "B1-A1",
    "A2-0", "B2-0", "A2-B2", "B2-A2",
    "A3-0", "B3-0", "A3-B3", "B3-A3",
    "A4-0", "B4-0", "A4-B4", "B4-A4",
    "A5-0", "B5-0", "A5-B5", "B5-A5",
    "A6-0", "B6-0", "A6-B6", "B6-A6",
    "Am1-B1", "Bm1-A1",
    "Am4-B1", "Bm4-A1",
    "Am6-B1", "Bm6-A1",
    "Am7-B1", "Bm7-A1"
)

all_alq_treatments <- c(
    "extraction_blank", "A1",
    "A2", "A3", "B1",
    "B2", "B3", "SW",
    "A4", "A5", "A6",
    "B4", "B5", "B6",
    "dropout_community"
)

pilot_comm_A <- c("ESL0825", "ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL0295", "ESL0185", "ESL0183", "ESL0184", "ESL0186")
pilot_comm_B <- c("ESL0820", "ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0197", "ESL0394", "ESL0263", "ESL0262", "ESL0350", "ESL0261")
comm_A1 <- c("ESL0825", "ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL0295", "ESL1028", "ESL0185", "ESL0183", "ESL0184", "ESL0186" )
comm_B1 <- c("ESL0820", "ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0197", "ESL0394", "ESL1028", "ESL0263", "ESL0835", "ESL0350", "ESL0261")
comm_A2 <- c("ESL0820", "ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0197", "ESL0295", "ESL1028", "ESL0185", "ESL0183", "ESL0184", "ESL0186" )
comm_B2 <- c("ESL0825", "ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL0394", "ESL1028", "ESL0263", "ESL0835", "ESL0350", "ESL0261")
comm_A3 <- c("ESL0820", "ESL0170", "ESL0824", "ESL0827", "ESL0819", "ESL0197", "ESL0394", "ESL1028", "ESL0185", "ESL0183", "ESL0350", "ESL0186" )
comm_B3 <- c("ESL0825", "ESL0822", "ESL0198", "ESL0199", "ESL0200", "ESL0197", "ESL0295", "ESL1028", "ESL0263", "ESL0835", "ESL0184", "ESL0261")
comm_A4 <- c("ESL0820", "ESL0822", "ESL0198", "ESL0827", "ESL0819", "ESL0197", "ESL0295", "ESL1028", "ESL0263", "ESL0183", "ESL0350", "ESL0186" )
comm_B4 <- c("ESL0825", "ESL0170", "ESL0824", "ESL0199", "ESL0200", "ESL0197", "ESL0394", "ESL1028", "ESL0185", "ESL0835", "ESL0184", "ESL0261")
comm_A5 <- c("ESL0820", "ESL0822", "ESL0198", "ESL0199", "ESL0200", "ESL0197", "ESL0295", "ESL1028", "ESL0263", "ESL0835", "ESL0184", "ESL0186" )
comm_B5 <- c("ESL0825", "ESL0170", "ESL0824", "ESL0827", "ESL0819", "ESL0197", "ESL0394", "ESL1028", "ESL0185", "ESL0183", "ESL0350", "ESL0261")
comm_A6 <- c("ESL0825", "ESL0170", "ESL0824", "ESL0199", "ESL0819", "ESL0197", "ESL0295", "ESL1028", "ESL0185", "ESL0835", "ESL0184", "ESL0261")
comm_B6 <- c("ESL0820", "ESL0822", "ESL0198", "ESL0827", "ESL0200", "ESL0197", "ESL0394", "ESL1028", "ESL0263", "ESL0183", "ESL0350", "ESL0186" )
comm_A7 <- c("ESL0825", "ESL0822", "ESL0824", "ESL0199", "ESL0819", "ESL0197", "ESL0394", "ESL1028", "ESL0185", "ESL0835", "ESL0184", "ESL0261")
comm_B7 <- c("ESL0820", "ESL0170", "ESL0198", "ESL0827", "ESL0200", "ESL0197", "ESL0295", "ESL1028", "ESL0263", "ESL0183", "ESL0350", "ESL0186" )
comm_Am1 <- c("ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL0295", "ESL1028", "ESL0185", "ESL0183", "ESL0184", "ESL0186" )
comm_Bm1 <- c("ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0197", "ESL0394", "ESL1028", "ESL0263", "ESL0835", "ESL0350", "ESL0261")
comm_Am4 <- c("ESL0825", "ESL0822", "ESL0824", "ESL0200", "ESL0197", "ESL0295", "ESL1028", "ESL0185", "ESL0183", "ESL0184", "ESL0186" )
comm_Bm4 <- c("ESL0820", "ESL0170", "ESL0198", "ESL0819", "ESL0197", "ESL0394", "ESL1028", "ESL0263", "ESL0835", "ESL0350", "ESL0261")
comm_Am6 <- c("ESL0825", "ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL1028", "ESL0185", "ESL0183", "ESL0184", "ESL0186" )
comm_Bm6 <- c("ESL0820", "ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0197", "ESL1028", "ESL0263", "ESL0835", "ESL0350", "ESL0261")
comm_Am7 <- c("ESL0825", "ESL0822", "ESL0824", "ESL0827", "ESL0200", "ESL0197", "ESL0295", "ESL1028", "ESL0183", "ESL0184", "ESL0186" )
comm_Bm7 <- c("ESL0820", "ESL0170", "ESL0198", "ESL0199", "ESL0819", "ESL0197", "ESL0394", "ESL1028", "ESL0835", "ESL0350", "ESL0261")


get_first_second_only_arriver <- function(my_treatment, my_strain) {
  if (my_strain == "ESL0353_ESL0185") {
    my_strain = "ESL0185"
  }
  if (my_strain == "ESL0183_ESL0262") {
    my_strain = "ESL0183"
  }
  if (!(my_strain %in% unique(c(comm_A1, comm_B1, pilot_comm_A, pilot_comm_B)))) {
    return(NA)
  }
  if (my_strain %in% c("ESL0197", "ESL1028")) {
    return("both")
  }
  if (my_treatment %in% all_alq_treatments) {
    return("aliquot")
  }
  if (my_treatment %in% c("MD")) {
    return("MD")
  }
  if (my_treatment %in% c("A1-0", "A1-B1")) {
    if (my_strain %in% comm_A1) {
      if (my_treatment == "A1-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "A1-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("B1-0", "B1-A1")) {
    if (my_strain %in% comm_B1) {
      if (my_treatment == "B1-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "B1-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("A2-0", "A2-B2")) {
    if (my_strain %in% comm_A2) {
      if (my_treatment == "A2-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "A2-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("B2-0", "B2-A2")) {
    if (my_strain %in% comm_B2) {
      if (my_treatment == "B2-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "B2-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("A3-0", "A3-B3")) {
    if (my_strain %in% comm_A3) {
      if (my_treatment == "A3-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "A3-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("B3-0", "B3-A3")) {
    if (my_strain %in% comm_B3) {
      if (my_treatment == "B3-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "B3-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("A4-0", "A4-B4")) {
    if (my_strain %in% comm_A4) {
      if (my_treatment == "A4-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "A4-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("B4-0", "B4-A4")) {
    if (my_strain %in% comm_B4) {
      if (my_treatment == "B4-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "B4-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("A5-0", "A5-B5")) {
    if (my_strain %in% comm_A5) {
      if (my_treatment == "A5-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "A5-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("B5-0", "B5-A5")) {
    if (my_strain %in% comm_B5) {
      if (my_treatment == "B5-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "B5-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("A6-0", "A6-B6")) {
    if (my_strain %in% comm_A6) {
      if (my_treatment == "A6-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "A6-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("B6-0", "B6-A6")) {
    if (my_strain %in% comm_B6) {
      if (my_treatment == "B6-0") {
        return("only")
      } else {
        return("first")
      }
    } else {
      if (my_treatment == "B6-0") {
        return("unadded")
      } else {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("Am1-B1", "Am4-B1", "Am6-B1", "Am7-B1")) {
    if (my_strain %in% comm_B1) {
      return("second")
    } else {
      if( my_treatment == "Am1-B1") {
        if (my_strain %in% comm_Am1) {
          return("first")
        } else {
          return("dropout")
        }
      }
      if( my_treatment == "Am4-B1") {
        if (my_strain %in% comm_Am4) {
          return("first")
        } else {
          return("dropout")
        }
      }
      if( my_treatment == "Am6-B1") {
        if (my_strain %in% comm_Am6) {
          return("first")
        } else {
          return("dropout")
        }
      }
      if( my_treatment == "Am7-B1") {
        if (my_strain %in% comm_Am7) {
          return("first")
        } else {
          return("dropout")
        }
      }
    }
  }
  if (my_treatment %in% c("Bm1-A1", "Bm4-A1", "Bm6-A1", "Bm7-A1")) {
    if (my_strain %in% comm_A1) {
      return("second")
    } else {
      if (my_treatment == "Bm1-A1") {
        if (my_strain %in% comm_Bm1) {
          return("first")
        } else {
          return("dropout")
        }
      }
      if (my_treatment == "Bm4-A1") {
        if (my_strain %in% comm_Bm4) {
          return("first")
        } else {
          return("dropout")
        }
      }
      if (my_treatment == "Bm6-A1") {
        if (my_strain %in% comm_Bm6) {
          return("first")
        } else {
          return("dropout")
        }
      }
      if (my_treatment == "Bm7-A1") {
        if (my_strain %in% comm_Bm7) {
          return("first")
        } else {
          return("dropout")
        }
      }
    }
  }
  if (my_treatment %in% c("AB-3")) {
    if (my_strain %in% pilot_comm_A) {
      return("first")
    } else {
      if (my_strain %in% pilot_comm_B) {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("BA-3")) {
    if (my_strain %in% pilot_comm_B) {
      return("first")
    } else {
      if (my_strain %in% pilot_comm_A) {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("AB-1")) {
    if (my_strain %in% pilot_comm_A) {
      return("first")
    } else {
      if (my_strain %in% pilot_comm_B) {
        return("second")
      }
    }
  }
  if (my_treatment %in% c("BA-1")) {
    if (my_strain %in% pilot_comm_B) {
      return("first")
    } else {
      if (my_strain %in% pilot_comm_A) {
        return("second")
      }
    }
  }
  if (my_treatment == "AB-0") {
    return("both")
  }
  if (my_treatment %in% c("A-0", "A-3_14", "A-3")) {
    if (my_strain %in% pilot_comm_A) {
      return("only")
    } else {
      if (my_strain %in% pilot_comm_B) {
        return("unadded")
      }
    }
  }
  if (my_treatment %in% c("B-0", "B-3_14", "B-3")) {
    if (my_strain %in% pilot_comm_B) {
      return("only")
    } else {
      if (my_strain %in% pilot_comm_A) {
        return("unadded")
      }
    }
  }
  return(NA)
}

# instead of "first" and "second" use "1_first" and "2_second"
get_first_second_only_arriver_mod <- function(my_treatment, my_strain) {
  if (my_strain == "ESL0353_ESL0185") {
    my_strain = "ESL0185"
  }
  if (my_strain == "ESL0183_ESL0262") {
    my_strain = "ESL0183"
  }
  if (!(my_strain %in% unique(c(comm_A1, comm_B1, pilot_comm_A, pilot_comm_B)))) {
    return(NA)
  }
  if (my_strain %in% c("ESL0197", "ESL1028")) {
    return("both")
  }
  if (my_treatment %in% all_alq_treatments) {
    return("aliquot")
  }
  if (my_treatment %in% c("MD")) {
    return("MD")
  }
  if (my_treatment %in% c("A1-0", "A1-B1")) {
    if (my_strain %in% comm_A1) {
      if (my_treatment == "A1-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "A1-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("B1-0", "B1-A1")) {
    if (my_strain %in% comm_B1) {
      if (my_treatment == "B1-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "B1-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("A2-0", "A2-B2")) {
    if (my_strain %in% comm_A2) {
      if (my_treatment == "A2-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "A2-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("B2-0", "B2-A2")) {
    if (my_strain %in% comm_B2) {
      if (my_treatment == "B2-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "B2-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("A3-0", "A3-B3")) {
    if (my_strain %in% comm_A3) {
      if (my_treatment == "A3-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "A3-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("B3-0", "B3-A3")) {
    if (my_strain %in% comm_B3) {
      if (my_treatment == "B3-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "B3-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("A4-0", "A4-B4")) {
    if (my_strain %in% comm_A4) {
      if (my_treatment == "A4-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "A4-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("B4-0", "B4-A4")) {
    if (my_strain %in% comm_B4) {
      if (my_treatment == "B4-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "B4-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("A5-0", "A5-B5")) {
    if (my_strain %in% comm_A5) {
      if (my_treatment == "A5-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "A5-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("B5-0", "B5-A5")) {
    if (my_strain %in% comm_B5) {
      if (my_treatment == "B5-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "B5-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("A6-0", "A6-B6")) {
    if (my_strain %in% comm_A6) {
      if (my_treatment == "A6-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "A6-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("B6-0", "B6-A6")) {
    if (my_strain %in% comm_B6) {
      if (my_treatment == "B6-0") {
        return("0_only")
      } else {
        return("1_first")
      }
    } else {
      if (my_treatment == "B6-0") {
        return("unadded")
      } else {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("Am1-B1", "Am4-B1", "Am6-B1", "Am7-B1")) {
    if (my_strain %in% comm_B1) {
      return("2_second")
    } else {
      if( my_treatment == "Am1-B1") {
        if (my_strain %in% comm_Am1) {
          return("1_first")
        } else {
          return("dropout")
        }
      }
      if( my_treatment == "Am4-B1") {
        if (my_strain %in% comm_Am4) {
          return("1_first")
        } else {
          return("dropout")
        }
      }
      if( my_treatment == "Am6-B1") {
        if (my_strain %in% comm_Am6) {
          return("1_first")
        } else {
          return("dropout")
        }
      }
      if( my_treatment == "Am7-B1") {
        if (my_strain %in% comm_Am7) {
          return("1_first")
        } else {
          return("dropout")
        }
      }
    }
  }
  if (my_treatment %in% c("Bm1-A1", "Bm4-A1", "Bm6-A1", "Bm7-A1")) {
    if (my_strain %in% comm_A1) {
      return("2_second")
    } else {
      if (my_treatment == "Bm1-A1") {
        if (my_strain %in% comm_Bm1) {
          return("1_first")
        } else {
          return("dropout")
        }
      }
      if (my_treatment == "Bm4-A1") {
        if (my_strain %in% comm_Bm4) {
          return("1_first")
        } else {
          return("dropout")
        }
      }
      if (my_treatment == "Bm6-A1") {
        if (my_strain %in% comm_Bm6) {
          return("1_first")
        } else {
          return("dropout")
        }
      }
      if (my_treatment == "Bm7-A1") {
        if (my_strain %in% comm_Bm7) {
          return("1_first")
        } else {
          return("dropout")
        }
      }
    }
  }
  if (my_treatment %in% c("AB-3")) {
    if (my_strain %in% pilot_comm_A) {
      return("1_first")
    } else {
      if (my_strain %in% pilot_comm_B) {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("BA-3")) {
    if (my_strain %in% pilot_comm_B) {
      return("1_first")
    } else {
      if (my_strain %in% pilot_comm_A) {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("AB-1")) {
    if (my_strain %in% pilot_comm_A) {
      return("1_first")
    } else {
      if (my_strain %in% pilot_comm_B) {
        return("2_second")
      }
    }
  }
  if (my_treatment %in% c("BA-1")) {
    if (my_strain %in% pilot_comm_B) {
      return("1_first")
    } else {
      if (my_strain %in% pilot_comm_A) {
        return("2_second")
      }
    }
  }
  if (my_treatment == "AB-0") {
    return("both")
  }
  if (my_treatment %in% c("A-0", "A-3_14", "A-3")) {
    if (my_strain %in% pilot_comm_A) {
      return("0_only")
    } else {
      if (my_strain %in% pilot_comm_B) {
        return("unadded")
      }
    }
  }
  if (my_treatment %in% c("B-0", "B-3_14", "B-3")) {
    if (my_strain %in% pilot_comm_B) {
      return("0_only")
    } else {
      if (my_strain %in% pilot_comm_A) {
        return("unadded")
      }
    }
  }
  return(NA)
}


pcoa_plot <- function(df_pcoa, metadata, variable, color_add=F, color_list, colname_in_metadata = "ID", shape_add = F, shape_var) {
          matrix <- as.matrix(df_pcoa)
          dist <- as.dist(matrix)
          res_pcoa <- pcoa(dist)
          ev1 <- res_pcoa$vectors[,1]
          ev2 <- res_pcoa$vectors[,2]
          df_pcoa_new <- data.frame(cbind(ev1,ev2))
          df_pcoa_new$Sample <- rownames(df_pcoa_new)
          rownames(df_pcoa_new) <- NULL
          df_pcoa_new <- left_join(df_pcoa_new, metadata, by = c("Sample" = colname_in_metadata))
          perc_axis <- round(((res_pcoa$values$Relative_eig[c(1,2)])*100), digits=1)
          axis_x_title <- paste0("PCo1 (",perc_axis[1],"%)")
          axis_y_title <- paste0("PCo2 (",perc_axis[2],"%)")
          if(color_add & shape_add) {
            p <- ggplot(df_pcoa_new, aes(x = ev1,
                                       y = ev2,
                                       shape = get(shape_var),
                                       colour = get(variable)))+
                geom_point(stat="identity", size=4) +
                  labs(x=axis_x_title, y = axis_y_title, color = variable, shape = shape_var) +
                    make_theme(setFill = F, setCol = F, guide_nrow = 7, leg_size = 10 ) +
                      scale_color_manual(values=color_list)
          } else {
            if (color_add) {
              p <- ggplot(df_pcoa_new, aes(x = ev1,
                                       y = ev2,
                                       colour = get(variable)))+
                geom_point(stat="identity", size=4, shape=19) +
                  labs(x=axis_x_title, y = axis_y_title, color = variable) +
                    make_theme(setFill = F, setCol = F, guide_nrow = 7, leg_size = 10 ) +
                      scale_color_manual(values=color_list)
            } else {
              p <- ggplot(df_pcoa_new, aes(x = ev1,
                                       y = ev2,
                                       color = get(variable)))+
                geom_point(stat="identity", size=4, shape=19) +
                  labs(x=axis_x_title, y = axis_y_title, color = variable) +
                    make_theme( max_colors = length(unique(df_pcoa_new[, variable])), guide_nrow = 7, leg_size = 10 ) 
            }
          }
          return(p)
}


# make treatment colors such that A1-0 and B1-0 are contrasting but A1-B1 and B1-A1 are similar to A1 and B1 each
# A1, A2, A3... all have the same color

TreatmentColors <- c(
  "A" = "#1f78b4",
  "B" = "#e31a1c",
  "AB" = "#a6cee3",
  "BA" = "#fb9a99",
  "AmB" = "#08519c",
  "BmA" = "#a50f15"
)

simplify_treatment <- function(my_treatment) {
  if (my_treatment %in% c("A1-0", "A2-0", "A3-0", "A4-0", "A5-0", "A6-0")) {
    return("A")
  }
  if (my_treatment %in% c("B1-0", "B2-0", "B3-0", "B4-0", "B5-0", "B6-0")) {
    return("B")
  }
  if (my_treatment %in% c("A1-B1", "A2-B2", "A3-B3", "A4-B4", "A5-B5", "A6-B6")) {
    return("AB")
  }
  if (my_treatment %in% c("B1-A1", "B2-A2", "B3-A3", "B4-A4", "B5-A5", "B6-A6")) {
    return("BA")
  }
  if (my_treatment %in% c("Am1-B1", "Am4-B1", "Am6-B1", "Am7-B1")) {
    return("AmB")
  }
  if (my_treatment %in% c("Bm1-A1", "Bm4-A1", "Bm6-A1", "Bm7-A1")) {
    return("BmA")
  }
  return(my_treatment)
}


# This function is from https://github.com/Russel88/MicEco/blob/91f8e6f5d67e0bfd018dc3b46da3994b1eadeb46/R/adonis_OmegaSq.R#L10
#' Calculate (partial) Omega-squared (effect-size calculation) for PERMANOVA and add it to the input object
#'
#' @param adonisOutput An adonis object
#' @param partial Should partial omega-squared be calculated (sample size adjusted). Default TRUE
#' @return Original adonis object with the (partial) Omega-squared values added
#' @import vegan
#' @export
adonis_OmegaSq <- function(adonisOutput, partial = TRUE){
    if(!(is(adonisOutput, "adonis") || is(adonisOutput, "anova.cca")))
        stop("Input should be an adonis object")
    if (is(adonisOutput, "anova.cca")) {
        aov_tab <- adonisOutput
        aov_tab$MeanSqs <- aov_tab$SumOfSqs / aov_tab$Df
        aov_tab$MeanSqs[length(aov_tab$Df)] <- NA
    } else {
        aov_tab <- adonisOutput$aov.tab
    }
    heading <- attr(aov_tab, "heading")
    MS_res <- aov_tab[pmatch("Residual", rownames(aov_tab)), "MeanSqs"]
    SS_tot <- aov_tab[rownames(aov_tab) == "Total", "SumsOfSqs"]
    N <- aov_tab[rownames(aov_tab) == "Total", "Df"] + 1
    if(partial){
        omega <- apply(aov_tab, 1, function(x) (x["Df"]*(x["MeanSqs"]-MS_res))/(x["Df"]*x["MeanSqs"]+(N-x["Df"])*MS_res))
        aov_tab$parOmegaSq <- c(omega[1:(length(omega)-2)], NA, NA)
    } else {
        omega <- apply(aov_tab, 1, function(x) (x["SumsOfSqs"]-x["Df"]*MS_res)/(SS_tot+MS_res))
        aov_tab$OmegaSq <- c(omega[1:(length(omega)-2)], NA, NA)
    }
    if (is(adonisOutput, "adonis"))
        cn_order <- c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2",
                      if (partial) "parOmegaSq" else "OmegaSq", "Pr(>F)")
    else
        cn_order <- c("Df", "SumOfSqs", "F", if (partial) "parOmegaSq" else "OmegaSq",
                      "Pr(>F)")
    aov_tab <- aov_tab[, cn_order]
    attr(aov_tab, "names") <- cn_order
    attr(aov_tab, "heading") <- heading
    if (is(adonisOutput, "adonis"))
        adonisOutput$aov.tab <- aov_tab
    else
        adonisOutput <- aov_tab
    return(adonisOutput)
}


strainShapesBySpec <- c(# Bifidobacterium
    "ESL0825" = 21, # Bifidobacterium apousia
    "ESL0820" = 24, # Bifidobacterium apousia
    "ESL0822" = 21, # Bifidobacterium asteroides
    "ESL0170" = 24, # Bifidobacterium asteroides
    "ESL0824" = 21, # Bifidobacterium polysaccharolyticum
    "ESL0198" = 24, # Bifidobacterium polysaccharolyticum
    "ESL0827" = 21, # Bifidobacterium sp1.
    "ESL0199" = 24, # Bifidobacterium sp1.
    "ESL0200" = 21, # Bifidobacterium sp2.
    "ESL0819" = 24, # Bifidobacterium sp2.
    "ESL0197" = 21, # Bifidobacterium coryneforme

    # Bombilactobacillus
    "ESL0295" = 21, # Bombilactobacillus mellis
    "ESL1028" = 21, # Bombilactobacillus mellifer

    # Lactobacillus
    "ESL0185" = 21, # Lactobacillus apis
    "ESL0353_ESL0185" = 21, # Lactobacillus apis
    "ESL0263" = 24, # Lactobacillus apis
    "ESL0183" = 21, # Lactobacillus helsingborgensis
    "ESL0183_ESL0262" = 21, # Lactobacillus helsingborgensis
    "ESL0835" = 24, # Lactobacillus helsingborgensis
    "ESL0184" = 21, # Lactobacillus melliventris
    "ESL0350" = 24, # Lactobacillus melliventris
    
    "ESL0394" = 23, # Lactobacillus melliventris
    
    "ESL0186" = 21, # Lactobacillus kullabergensis
    "ESL0261" = 24  # Lactobacillus kullabergensis
)

HostColors <- c(
  "Apis mellifera" = "#1f78b4",
  "Apis cerana" = "#e31a1c",
  "Apis dorsata" = "#6a3d9a",
  "Apis florea" = "#33a02c",
  "Apis andreniformis" = "#ff7f00"
)
