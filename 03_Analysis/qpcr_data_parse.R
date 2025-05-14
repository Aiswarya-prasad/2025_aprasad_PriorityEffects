source('03_Analysis/Utilities.R', chdir = TRUE)

# metadata_path <- paste0(prefix_dir, "/spirit/D2c/aprasad/20220921_aprasad_PriorityEffectsExperimentPilot/03_PilotExperiment/samples_metadata.xlsx")
# df_meta <- read_excel(metadata_path, sheet = 1)


prefix_dir <- "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/03_Analysis/PriorityEffectsExperiment/RawData/SampleqPCR/"

qpcr_results <- c(
  paste0(prefix_dir, "plate1-1.xls"),
  paste0(prefix_dir, "plate1-2.xls"),
  paste0(prefix_dir, "plate1-3.xls"),
  paste0(prefix_dir, "plate2-1.xls"),
  paste0(prefix_dir, "plate2-2.xls"),
  paste0(prefix_dir, "plate3-1.xls"),
  paste0(prefix_dir, "plate3-2.xls"),
  paste0(prefix_dir, "plate3-3.xls"),
  paste0(prefix_dir, "plate4-1.xls"),
  paste0(prefix_dir, "plate4-2.xls"),
  paste0(prefix_dir, "plate5-1.xls"),
  paste0(prefix_dir, "plate5-2.xls"),
  paste0(prefix_dir, "plate5-3.xls"),
  paste0(prefix_dir, "plate6-1.xls"),
  paste0(prefix_dir, "plate6-2.xls"),
  paste0(prefix_dir, "plate6-3.xls")
)

qpcr_df <- data.frame()
plate_num <- 1
for(x in qpcr_results){
  # read just the rows and columns with the values and leave out metadata cells
  # Limits the cells from rows=(36 to 132) and columns=(2 to 11)
  t <- read_excel(x, sheet = 1, range = cell_limits(c(36,2), c(132,11)))
  t <- t[,c(1,3,4,5,8,9,10)]
  names(t)[1:5] <- c("Well", "Sample.Name","Target.Name","Task","Ct")
  t$Ct <- as.numeric(t$Ct) # NAs will be introduced at the place of 'undetermined' values
  t$Target.Name <- factor(t$Target.Name)
  t$Sample.Name <- factor(t$Sample.Name)
  # t$sd <- factor(t$sd)
  t <- as.data.frame(t)
  # add in the plate number
  t <- cbind(t, "Plate" = as.factor(rep(plate_num, dim(t)[1])))
  # bind the previous plate values with the next before moving to the next file
  qpcr_df <- rbind(qpcr_df, t)
  plate_num = plate_num + 1
}

qpcr_df_info <- qpcr_df %>%
          mutate(Target.Name = ifelse(Target.Name == "UV_0356/0772", "Bacteria", "Host")) %>%
            mutate(Sample.Name = as.character(Sample.Name)) %>%
                mutate(ID = ifelse(Task == "NTC", "neg", Sample.Name)) %>%
                    select(!Task)
                # mutate(ID = as.character(Sample.Name))
colnames(qpcr_df_info) <- c("Well", "Sample.Name", "Target.Name", "Ct", "Ct_Mean", "Ct_SD", "Plate", "ID")

control_values <- qpcr_df_info %>%
    filter(grepl("neg|pos", ID)) %>%
    mutate(ID = ifelse(grepl("positive", ID), "pos", ID)) %>%
    mutate(ID = ifelse(grepl("negative", ID), "neg", ID)) %>%
        select(ID, Target.Name, Ct, Ct_Mean, Plate) %>%
            group_by(Target.Name, ID, Plate) %>%
            mutate(Mean_control = mean(Ct, na.rm = T)) %>%
            select(Plate, Target.Name, ID, Mean_control) %>%
            unique() %>%
            mutate(Type = paste0(Target.Name, "_", ID)) %>%
            ungroup() %>%
            select(Plate, Mean_control, Type)

qpcr_df_bac <- qpcr_df_info %>%
                filter(!grepl("Control", ID)) %>%
                filter(!grepl("neg|pos", ID)) %>%
                    filter(Target.Name == "Bacteria") %>%
                        mutate(Ct_Mean_bac = Ct_Mean) %>%
                        select(!Ct_Mean) %>%
                        mutate(Ct_SD_bac = Ct_SD) %>%
                        select(!Ct_SD) %>%
                        select(!c(Well, Sample.Name, Target.Name, Ct)) %>%
                        unique()

############
# using values from NAS std curve 
# smb://nas.unil.ch/DMF/GROUPS/gr_Engel/lab_resources/Protocols/qPCR/qPCR_analysis_workflow/StdCurves/StdCurves.xlsx copied to
# done in Berta on 12-11-2022
# All the qPCRs were done on Berta and using univ primers (0356/0772)
############
UV_eff <- 1.971
UV_int <- 36.94109669
Actin_eff <- 2.001
Actin_int <- 36.71918987
elution_vol <- 30
dilution_factor <- 10
# DNA extraction protocol uses 80 out of 360 uL gut homogenate
# so copies obtained need to be multiplied by 360/80 = 4.5
extraction_frac <- 360/80
LOD <- 1000


############
df_qpcr_counts <- qpcr_df_bac %>%
                    filter(!is.na(ID)) %>%
                    mutate(bac_copies_raw = UV_eff^(UV_int - Ct_Mean_bac)) %>%
                    # copies per gut accounts for fraction of the homogenate, volume in which it was elutes times dilution factor
                    # because if the gut was eluted in 30 ul and 1 ul was used for qPCR in a 1:10 dilution, then the copies in the gut
                    # would be 30 times the copies in the qPCR what was in 1 uL used for qPCR
                    mutate(bac_copies = bac_copies_raw*elution_vol*dilution_factor*extraction_frac)

system("mkdir -p 03_Analysis/PriorityEffectsExperiment/01_qpcrdata")
write.csv(df_qpcr_counts, "03_Analysis/PriorityEffectsExperiment/01_qpcrdata/qPCR_counts.csv", row.names = F, quote = F)
