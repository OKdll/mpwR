
#import
import_mpwR <- function(input_file, ...) {

  data.table::fread(input_file, ...) %>%
    tibble::as_tibble()

}

# Prepare the input data
# Input data will be renamed and default filtering will be applied
# Default filtering- MaxQuant: no potential contaminants and no reverse sequences; DIA-NN: none; Spectronaut: EG.Identified only TRUE; PD: Only High Confidence identifications

prepare_input <- function(input_df,
                          software = c("MaxQuant", "DIA-NN", "Spectronaut", "PD"),
                          MaxQuant_addon = c("evidence", "peptide", "proteingroup"),
                          PD_addon = c("psm", "peptide", "protein", "proteingroup")) {

  #handle global vars
  . <- NULL

  #**dependencies
  if (software != "MaxQuant" & software != "DIA-NN" & software != "Spectronaut" & software != "PD") {
    stop("Did you spell the software input wrong? Choose between: MaxQuant, DIA-NN, Spectronaut, PD")
  }
  ##

  #****************#
  #*Add all columns for identifying file-software match and columns for downstream analysis:

  cols_MQ_ev <- c("Raw file", "Proteins", "Modified sequence", "Sequence", "Missed cleavages", "Charge", "Retention time", "Potential contaminant", "Reverse")
  cols_MQ_pep <- c("Sequence", "Missed cleavages", "Last amino acid", "Amino acid after", "Potential contaminant", "Reverse") #"Last amino acid" and "amino acid after" to clearly identify data as peptide.txt - currently not used for downstream analysis #requires intensity columns for tidying #requires LFQ intensity for CV LFQ calc
  cols_MQ_pg <- c("Protein IDs", "Majority protein IDs", "Peptide counts (all)", "Potential contaminant", "Reverse") #"Majority protein IDs" and "Peptide counts (all)" to clearly identify peptideGroups.txt - currently not used for downstream analysis #requires intensity columns for tidying #requires LFQ intensity for CV LFQ calc

  cols_DIANN <- c("Run", "Protein.Group", "Protein.Ids", "Modified.Sequence", "Stripped.Sequence", "Precursor.Id", "Precursor.Charge", "RT", "PG.MaxLFQ")

  cols_spec <- c("R.FileName", "PG.ProteinGroups", "PEP.StrippedSequence", "EG.ModifiedPeptide", "EG.PrecursorId", "EG.Identified", "EG.ApexRT", "FG.Charge", "PG.Quantity", "PEP.Quantity", "PEP.NrOfMissedCleavages") # currently not used: PEP.IsProteotypic, EG.IsDecoy

  cols_PD_psm <- c("Confidence", "Spectrum File", "Protein Accessions", "Annotated Sequence", "Modifications", "Number of Missed Cleavages", "Charge", "RT in min") # currently not used Apex RT in min
  cols_PD_pep <- c("Peptide Groups Peptide Group ID", "Confidence", "Sequence", "Modifications", "Number of Missed Cleavages") # also contains(Found in Sample) #currently not used Annotated Sequence #"Peptide Groups Peptide Group" ID to clearly identify peptideGroups.txt - not used for downstream analysis
  cols_PD_prot <- c("Proteins Unique Sequence ID", "Accession", "Description") # also contains(Found in Sample) #"Proteins Unique Sequence ID" and "Description" to cleary identify Proteins.txt - not used for downstream analysis
  cols_PD_pg <- c("Protein Groups Protein Group ID", "Group Description", "Number of Proteins", "Number of Unique Peptides") # also contains(Found in Sample) #"Protein Groups Protein Group ID","Number of Proteins", "Number of Unique Peptides" to clearly identify ProteinGroups.txt - not used for downstream analysis
  #***************

  if (software == "MaxQuant") {

    #dependencies
    if (MaxQuant_addon != "evidence" & MaxQuant_addon != "peptide" & MaxQuant_addon != "proteingroup") {
      stop("Did you spell the input for MaxQuant_addon wrong? Choose between: evidence, peptide, proteingroup")
    }
    #**

    if (MaxQuant_addon == "evidence") {

      #Column check
      if (sum(cols_MQ_ev %in% colnames(input_df)) != length(cols_MQ_ev)) {
        message <- paste("For MaxQuant evidence.txt - the following column need to be present: ", cols_MQ_ev, "\n")
        message(message)
        stop("Not all required columns present in submitted data.")
      }
      #**

      output_df <- input_df %>%
        dplyr::filter(.data$`Potential contaminant` != "+", .data$Reverse != "+") %>%
        dplyr::mutate("Run_mpwR" = .data$`Raw file`,
                      "Protein.IDs_mpwR" = .data$Proteins,
                      "Peptide.IDs_mpwR" = .data$`Modified sequence`,
                      "Stripped.Sequence_mpwR" = .data$Sequence,
                      "Missed.Cleavage_mpwR" = .data$`Missed cleavages`,
                      "Retention.time_mpwR" = .data$`Retention time`,
                      "Precursor.Charge_mpwR" = .data$Charge) %>%
        dplyr::filter(.data$Peptide.IDs_mpwR != "") %>%
        tidyr::unite(., col = "Precursor.IDs_mpwR", c("Peptide.IDs_mpwR", "Precursor.Charge_mpwR"), sep = "", remove = FALSE)


    } else if (MaxQuant_addon == "peptide") {

      #Column check
      if (sum(cols_MQ_pep %in% colnames(input_df)) != length(cols_MQ_pep)) {
        message <- paste("For MaxQuant peptides.txt - the following column need to be present: ", cols_MQ_pep, "\n")
        stop("Not all required columns present in submitted data.")
      }
      #**

      output_df <- input_df %>%
        dplyr::filter(.data$`Potential contaminant` != "+", .data$Reverse != "+") %>%
        dplyr::mutate(
          "Stripped.Sequence_mpwR" = .data$Sequence,
          "Missed.Cleavage_mpwR" = .data$`Missed cleavages`
        )

    } else if (MaxQuant_addon == "proteingroup") {

      #Column check
      if (sum(cols_MQ_pg %in% colnames(input_df)) != length(cols_MQ_pg)) {
        message <- paste("For MaxQuant proteinGroups.txt input - the following column need to be present: ", cols_MQ_pg, "\n")
        message(message)
        stop("Not all required columns present in submitted data.")
      }
      #**

      output_df <- input_df %>%
        dplyr::filter(.data$`Potential contaminant` != "+", .data$Reverse != "+") %>%
        dplyr::mutate("ProteinGroup.IDs_mpwR" = .data$`Protein IDs`)

    }

  } else if (software == "DIA-NN") {

    #Column check
    if (sum(cols_DIANN %in% colnames(input_df)) != length(cols_DIANN)) {
      message <- paste("For DIA-NN input - the following column need to be present: ", cols_DIANN, "\n")
      message(message)
      stop("Not all required columns present in submitted data.")
    }
    #**

    output_df <- input_df %>%
      dplyr::mutate(
        "Run_mpwR" = .data$Run,
        "Stripped.Sequence_mpwR" = .data$Stripped.Sequence,
        "ProteinGroup.IDs_mpwR" = .data$Protein.Group,
        "Protein.IDs_mpwR" = .data$Protein.Ids,
        "Peptide.IDs_mpwR" = .data$Modified.Sequence,
        "Precursor.IDs_mpwR" = .data$Precursor.Id,
        "ProteinGroup_LFQ_mpwR" = .data$PG.MaxLFQ,
        "Precursor.Charge_mpwR" = .data$Precursor.Charge,
        "Retention.time_mpwR" = .data$RT)

  } else if (software == "Spectronaut") {

    #Column check
    if (sum(cols_spec %in% colnames(input_df)) != length(cols_spec)) {
      message <- paste("For Spectronaut input - the following column need to be present: ", cols_spec, "\n")
      message(message)
      stop("Not all required columns present in submitted data.")
    }
    #**

    output_df <- input_df %>%
      dplyr::filter(.data$EG.Identified == "TRUE") %>%
      dplyr::mutate("Run_mpwR" = .data$R.FileName,
                    "ProteinGroup.IDs_mpwR" = .data$PG.ProteinGroups,
                    "Peptide.IDs_mpwR" = .data$EG.ModifiedPeptide,
                    "Precursor.IDs_mpwR" = .data$EG.PrecursorId,
                    "Stripped.Sequence_mpwR" = .data$PEP.StrippedSequence,
                    "Precursor.Charge_mpwR" = .data$FG.Charge,
                    "Missed.Cleavage_mpwR" = .data$PEP.NrOfMissedCleavages,
                    "Retention.time_mpwR" = .data$EG.ApexRT,
                    "ProteinGroup_LFQ_mpwR" = .data$PG.Quantity, #not not necessarily LFQ quantity - depends on analysis settings in spectronaut
                    "Peptide_LFQ_mpwR" = .data$PEP.Quantity) #not necessarily LFQ quantity - depends on analysis settings in spectronaut

  } else if (software == "PD") {

    #dependencies
    if (PD_addon != "psm" & PD_addon != "protein" & PD_addon != "peptide" & PD_addon != "proteingroup") {
      stop("Did you spell the input for PD_addon wrong? Choose between: psm, evidence, peptide, proteingroup")
    }
    #**

    if (PD_addon == "psm") {

      #Column check
      if (sum(cols_PD_psm %in% colnames(input_df)) != length(cols_PD_psm)) {
        message <- paste("For PD PSMs.txt - the following column need to be present: ", cols_PD_psm, "\n")
        message(message)
        stop("Not all required columns present in submitted data.")
      }
      #**

      output_df <- input_df %>%
        dplyr::filter(.data$Confidence == "High") %>% #other: "Medium" - "Low"
        dplyr::mutate("Run_mpwR" = .data$`Spectrum File`,
                      "ProteinGroup.IDs_mpwR" = .data$`Protein Accessions`,
                      "Peptide.IDs_mpwR" = .data$`Annotated Sequence`,
                      "Details_Mods_mpwR" = .data$Modifications,
                      "Missed.Cleavage_mpwR" = .data$`Number of Missed Cleavages`,
                      "Precursor.Charge_mpwR" = .data$Charge,
                      "Retention.time_mpwR" = .data$`RT in min`) %>%
        tidyr::unite(., col = "Precursor.IDs_mpwR", c("Peptide.IDs_mpwR", "Precursor.Charge_mpwR"), sep = "", remove = FALSE)


    } else if (PD_addon == "peptide") {

      #Column check
      if (sum(cols_PD_pep %in% colnames(input_df)) != length(cols_PD_pep) | sum(stringr::str_detect(colnames(input_df), pattern = "Found in Sample")) == 0) {
        message <- paste("For PD PeptideGroups.txt - the following column need to be present: ", cols_PD_pep, "\n")
        message("For PD PeptideGroups.txt - the following column(s) need to be present: Found in Sample")
        message(message)
        stop("Not all required columns present in submitted data.")
      }
      #**

      output_df <- input_df %>%
        dplyr::filter(.data$Confidence == "High") %>% #other: "Medium" - "Low"
        dplyr::mutate(
          "Stripped.Sequence_mpwR" = .data$`Sequence`,
          "Details_Mods_mpwR" = .data$Modifications,
          "Missed.Cleavage_mpwR" = .data$`Number of Missed Cleavages`
        ) %>%
        tidyr::pivot_longer(cols = contains("Found in Sample"), names_to = "Run_mpwR", values_to = "Found.in.Sample_values_mpwR") %>% #tidy #Run needed for ID_Report - generate_level_count
        dplyr::filter(.data$Found.in.Sample_values_mpwR == "High") #"Not Found" - "Peak Found" - "Medium" - "Low"


    } else if (PD_addon == "protein") {

      #Column check
      if (sum(cols_PD_prot %in% colnames(input_df)) != length(cols_PD_prot) | sum(stringr::str_detect(colnames(input_df), pattern = "Found in Sample")) == 0) {
        message <- paste("For PD Proteins.txt - the following column need to be present: ", cols_PD_prot, "\n")
        message("For PD Proteins.txt - the following column(s) need to be present: Found in Sample")
        message(message)
        stop("Not all required columns present in submitted data.")
      }
      #**

      output_df <- input_df %>%
        dplyr::filter(.data$`Protein FDR Confidence Combined` == "High") %>% #other: "Medium" - "Low"
        dplyr::mutate(
          "Protein.IDs_mpwR" = .data$Accession) %>%
        tidyr::pivot_longer(cols = contains("Found in Sample"), names_to = "Run_mpwR", values_to = "Found.in.Sample_values_mpwR") %>% #tidy #Run needed for ID_Report - generate_level_count
        dplyr::filter(.data$Found.in.Sample_values_mpwR == "High") #"Not Found" - "Peak Found" - "Medium" - "Low"

    } else if (PD_addon == "proteingroup") {

      #Column check####
      if (sum(cols_PD_pg %in% colnames(input_df)) != length(cols_PD_pg) | sum(stringr::str_detect(colnames(input_df), pattern = "Found in Sample")) == 0) {
        message <- paste("For PD ProteinGroups.txt - the following column need to be present: ", cols_PD_pg, "\n")
        message("For PD ProteinGroups.txt - the following column(s) need to be present: Found in Sample")
        message(message)
        stop("Not all required columns present in submitted data.")
      }
      ##############

      output_df <- input_df %>%
        dplyr::rename(
          "ProteinGroup.IDs_mpwR" = .data$`Protein Groups Protein Group ID`) %>% #, #only number! - PD Manual page 306
        tidyr::pivot_longer(cols = contains("Found in Sample"), names_to = "Run_mpwR", values_to = "Found.in.Sample_values_mpwR") %>% #tidy #Run needed for ID_Report - generate_level_count
        dplyr::filter(.data$Found.in.Sample_values_mpwR == "High") #"Not Found" - "Peak Found" - "Medium" - "Low"
    }
  }
  return(output_df)
}
