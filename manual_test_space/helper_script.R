devtools::install_github("Thatguy027/goldengateR", force = TRUE)

install.packages(                                                                                                                                                          
  "/Users/Stefan/UCLA/Projects/synthetic_yeast/project_plan/GG_Rpackage/goldengateR",                                                                                      
  repos = NULL, type = "source"                                                                                                                                            
)    

library(goldengateR)


setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

set_plasmid_dir("/Users/Stefan/UCLA/Projects/synthetic_yeast/project_plan/GG_Rpackage/goldengateR/manual_test_space/xylose_parts/")

assemble_interactive(enzyme = "BsaI", output_dir = "./plasmids")  


entry_vec <- load_plasmid("pYTK033")
backbone  <- digest_and_select(entry_vec, "BsaI", select = "largest")

# 2. PCR-amplify your gene with BsaI primer tails, digest, take the insert
amp    <- read_part("pcr_products/xylb_pcr.fa", topology = "linear")
insert <- digest_and_select(amp, "BsaI", select = "largest")

# 3. Ligate
asm <- ligate(list(backbone, insert),
              enzyme = "BsaI",
              name   = "xylb_part",
              output = "xylose_parts/xylb_part.gb")

# 2. PCR-amplify your gene with BsaI primer tails, digest, take the insert
amp    <- read_part("pcr_products/gal2mut_pcr.fa", topology = "linear")
insert <- digest_and_select(amp, "BsaI", select = "largest")

# 3. Ligate
asm <- ligate(list(backbone, insert),
              enzyme = "BsaI",
              name   = "gal2_part",
              output = "xylose_parts/gal2_part.gb")


results <- assemble_from_table("assemblies_example.tsv", output_dir = "xylose_parts")


assemble_interactive(enzyme = "BsmBI", output_dir = "/Users/Stefan/UCLA/Projects/synthetic_yeast/project_plan/GG_Rpackage/goldengateR/manual_test_space/xylose_parts")  
