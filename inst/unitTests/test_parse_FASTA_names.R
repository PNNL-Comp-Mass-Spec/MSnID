library("MSnID")

test_parse_FASTA_names <- function() {
  path_to_FASTA <- system.file("extdata/uniprot_rat_small.fasta.gz",
                               package = "MSnID")
  out <- parse_FASTA_names(path_to_FASTA)
  
  # Reference data.frame
  ref <- data.frame(
    feature = c("sp|P61943|SIA10_RAT", "sp|P56576|UH11_RAT", 
                "sp|Q920B6-4|KCNK2_RAT", "sp|Q9WUW9|S1C2A_RAT", 
                "sp|Q6PCU8|NDUV3_RAT", "sp|O35509|RB11B_RAT"),
    database = "sp", 
    unique_id = c("P61943", "P56576", "Q920B6-4", "Q9WUW9", "Q6PCU8", "O35509"),
    uniprot_acc = c("P61943", "P56576", "Q920B6", "Q9WUW9", "Q6PCU8", "O35509"),
    isoform = c(rep(NA, 2), 4, rep(NA, 3)),
    entry_name = c("SIA10_RAT", "UH11_RAT", "KCNK2_RAT", "S1C2A_RAT", 
                   "NDUV3_RAT", "RB11B_RAT"),
    description = c(
      "Type 2 lactosamine alpha-2,3-sialyltransferase", 
      "Unknown protein from spot P11 of 2D-PAGE of heart tissue (Fragment)",
      "Isoform 4 of Potassium channel subfamily K member 2", 
      "Sulfotransferase 1C2A",
      "NADH dehydrogenase [ubiquinone] flavoprotein 3, mitochondrial", 
      "Ras-related protein Rab-11B"
    ),
    organism = "Rattus norvegicus",
    organism_id = 10116,
    gene = c("St3gal6", NA, "Kcnk2", "Sult1c2a", "Ndufv3", "Rab11b"),
    protein_existence = c(2, 1, NA, 1, 3, 1),
    sequence_version = c(1, 1, NA, 2, 1, 4)
  )
  
  # Subset and check
  set.seed(100)
  idx <- sample(1:100, 6)
  out_sub <- out[idx, ]
  rownames(out_sub) <- NULL
  
  checkEquals(out_sub, ref)
}

