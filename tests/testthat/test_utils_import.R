
test_that("prepare_input works", {

#MaxQuant
   #evidence
   data <- tibble::tibble(
      "Raw file" = c("R01", "R02", "R03"),
      "Proteins" = c("A", "B", "B"),
      "Modified sequence" = c("AMod", "BMod", "BMod"),
      "Sequence" = c("A", "B", "B"),
      "Missed cleavages" = c(1, 2, 2),
      "Charge" = c(2, 3, 3),
      "Retention time" = c(3, 5, 5.2),
      "Potential contaminant" = c("", "+", ""),
      "Reverse" = c("", "", "+")
   )

  output <- prepare_input(input_df = data, software = "MaxQuant", MaxQuant_addon = "evidence")

  expect_s3_class(output, "tbl")
  expect_equal(nrow(output), 1)
  expect_equal(ncol(output), 17)
  expect_equal(output$`Modified sequence`, "AMod")
  expect_equal(output$`Raw file`, "R01")
  expect_error(prepare_input(input_df = data[, -1], software = "MaxQuant", MaxQuant_addon = "evidence"), "Not all required columns present in submitted data.")

  #peptide
  data <- tibble::tibble(
     "Sequence" = c("A", "B", "B"),
     "Missed cleavages" = c(1, 2, 2),
     "Last amino acid" = c("A", "B", "B"),
     "Amino acid after" = c("A", "B", "B"),
     "Potential contaminant" = c("", "+", ""),
     "Reverse" = c("", "", "+")
  )

  output <- prepare_input(input_df = data, software = "MaxQuant", MaxQuant_addon = "peptide")

  expect_s3_class(output, "tbl")
  expect_equal(nrow(output), 1)
  expect_equal(ncol(output), 8)
  expect_equal(output$Sequence, "A")
  expect_equal(output$Stripped.Sequence_mpwR, "A")
  expect_error(prepare_input(input_df = data[, -1], software = "MaxQuant", MaxQuant_addon = "peptide"), "Not all required columns present in submitted data.")

  #proteingroup
  data <- tibble::tibble(
     "Protein IDs" = c("A", "B", "B"),
     "Majority protein IDs" = c("A", "B", "B"),
     "Peptide counts (all)" = c(1, 2, 2),
     "Potential contaminant" = c("", "+", ""),
     "Reverse" = c("", "", "+")
  )

  output <- prepare_input(input_df = data, software = "MaxQuant", MaxQuant_addon = "proteingroup")

  expect_s3_class(output, "tbl")
  expect_equal(nrow(output), 1)
  expect_equal(ncol(output), 6)
  expect_equal(output$ProteinGroup.IDs_mpwR, "A")
  expect_equal(output$`Protein IDs`, "A")
  expect_error(prepare_input(input_df = data[, -1], software = "MaxQuant", MaxQuant_addon = "proteingroup"), "Not all required columns present in submitted data.")

#DIA-NN
  data <- tibble::tibble(
    Protein.Group = c("A", "B", "B"),
    Precursor.Id = c("AAAATGTIFTFR2", "AEDTAVYYC(UniMod:4)AK2", "AAC(Dummy_Modification)LLPK2"),
    Run = c("R01", "R02", "R03"),
    Stripped.Sequence = c("AAAATGTIFTFR", "AEDTAVYYCAK", "AACLLPK"),
    Protein.Ids = c("A", "B", "B"),
    Modified.Sequence = c("AAAATGTIFTFR2", "AEDTAVYYC(UniMod:4)AK2", "AAC(Dummy_Modification)LLPK2"),
    PG.MaxLFQ = c(5, 5, 7),
    Precursor.Charge = c(2, 2, 2),
    RT = c(2, 3, 3.6)
  )

  output <- prepare_input(input_df = data, software = "DIA-NN")

  expect_s3_class(output, "tbl")
  expect_equal(nrow(output), 3)
  expect_equal(ncol(output), 18)
  expect_equal(output$ProteinGroup.IDs_mpwR, c("A", "B", "B"))
  expect_equal(output$Protein.Ids, c("A", "B", "B"))
  expect_error(prepare_input(input_df = data[, -1], software = "DIA-NN"), "Not all required columns present in submitted data.")

#Spectronaut
  data <- tibble::tibble(
    EG.Identified = c(TRUE, TRUE, FALSE),
    R.FileName = c("R01", "R02", "R03"),
    PG.ProteinGroups = c("A", "B", "B"),
    EG.ModifiedPeptide = c("AMod", "BMod", "BMod"),
    EG.PrecursorId = c("AMod2", "BMod2", "BMod2"),
    PEP.StrippedSequence = c("A", "B", "B"),
    FG.Charge = c(2, 2, 2),
    PEP.NrOfMissedCleavages = c(0, 1, 1),
    EG.ApexRT = c(3, 3, 4),
    PG.Quantity = c(5, 5, 7),
    PEP.Quantity = c(6, 6, 8)
  )

  output <- prepare_input(input_df = data, software = "Spectronaut")

  expect_s3_class(output, "tbl")
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), 21)
  expect_equal(output$ProteinGroup.IDs_mpwR, c("A", "B"))
  expect_equal(output$PG.ProteinGroups, c("A", "B"))
  expect_error(prepare_input(input_df = data[, -1], software = "DIA-NN"), "Not all required columns present in submitted data.")

#PD
 #psm
 data <- tibble::tibble(
   Confidence = c("High", "High", "Low"),
   `Spectrum File` = c("R01", "R02", "R03"),
   `Protein Accessions` = c("A", "B", "B"),
   `Annotated Sequence` = c("AMod", "BMod", "BMod"),
    Modifications = c("Mod", "Mod", "Mod"),
    `Number of Missed Cleavages` = c(2, 1, 0),
    Charge = c(2, 2, 2),
    `RT in min` = c(3, 3, 3)
 )

 output <- prepare_input(input_df = data, software = "PD", PD_addon = "psm")

 expect_s3_class(output, "tbl")
 expect_equal(nrow(output), 2)
 expect_equal(ncol(output), 16)
 expect_equal(output$ProteinGroup.IDs_mpwR, c("A", "B"))
 expect_equal(output$Precursor.IDs_mpwR, c("AMod2", "BMod2"))
 expect_error(prepare_input(input_df = data[, -1], software = "PD", PD_addon = "psm"), "Not all required columns present in submitted data.")

 #petpide
 data <- tibble::tibble(
   `Peptide Groups Peptide Group ID` = c(1, 2, 2),
   Confidence = c("High", "High", "Low"),
   `Sequence` = c("A", "B", "B"),
   Modifications = c("AMod", "BMod", "BMod"),
   `Number of Missed Cleavages` = c(0, 0, 1),
   `Found in Sample R01` = c("High", "High", "Low"),
   `Found in Sample R02` = c("High", "Low", "Low")
 )

 data2 <- data %>% select(-contains("Found"))

 output <- prepare_input(input_df = data, software = "PD", PD_addon = "peptide")

 expect_s3_class(output, "tbl")
 expect_equal(nrow(output), 3)
 expect_equal(ncol(output), 10)
 expect_equal(output$Run_mpwR, c("Found in Sample R01", "Found in Sample R02", "Found in Sample R01"))
 expect_equal(output$Stripped.Sequence_mpwR, c("A", "A", "B"))
 expect_error(prepare_input(input_df = data[, -2], software = "PD", PD_addon = "peptide"), "Not all required columns present in submitted data.")
 expect_error(prepare_input(input_df = data2, software = "PD", PD_addon = "peptide"), "Not all required columns present in submitted data.")

 #protein
 data <- tibble::tibble(
   `Proteins Unique Sequence ID` = c(1, 2, 2),
    Description = c("A", "B", "B"),
   `Protein FDR Confidence Combined` = c("High", "High", "Low"),
    Accession = c("A", "B", "B"),
   `Found in Sample R01` = c("High", "High", "Low"),
   `Found in Sample R02` = c("High", "Low", "Low")
 )

 data2 <- data %>% select(-contains("Found"))

 output <- prepare_input(input_df = data, software = "PD", PD_addon = "protein")

 expect_s3_class(output, "tbl")
 expect_equal(nrow(output), 3)
 expect_equal(ncol(output), 7)
 expect_equal(output$Run_mpwR, c("Found in Sample R01", "Found in Sample R02", "Found in Sample R01"))
 expect_equal(output$Protein.IDs_mpwR, c("A", "A", "B"))
 expect_error(prepare_input(input_df = data[, -2], software = "PD", PD_addon = "protein"), "Not all required columns present in submitted data.")
 expect_error(prepare_input(input_df = data2, software = "PD", PD_addon = "protein"), "Not all required columns present in submitted data.")

 #proteingroup
 data <- tibble::tibble(
   `Group Description` = c("A", "B", "B"),
   `Number of Proteins` = c(2, 2, 3),
   `Number of Unique Peptides` = c(3, 1, 4),
   `Protein Groups Protein Group ID` = c(1, 2, 2),
   `Found in Sample R01` = c("High", "High", "Low"),
   `Found in Sample R02` = c("High", "Low", "Low")
 )

 data2 <- data %>% select(-contains("Found"))

 output <- prepare_input(input_df = data, software = "PD", PD_addon = "proteingroup")

 expect_s3_class(output, "tbl")
 expect_equal(nrow(output), 3)
 expect_equal(ncol(output), 6)
 expect_equal(output$Run_mpwR, c("Found in Sample R01", "Found in Sample R02", "Found in Sample R01"))
 expect_equal(output$ProteinGroup.IDs_mpwR, c(1, 1, 2))
 expect_error(prepare_input(input_df = data[, -2], software = "PD", PD_addon = "proteingroup"), "Not all required columns present in submitted data.")
 expect_error(prepare_input(input_df = data2, software = "PD", PD_addon = "proteingroup"), "Not all required columns present in submitted data.")

})

