# Workflow test: build a Type 3 entry plasmid by gel-purifying the
# CamR+ColE1 backbone of an existing entry vector and ligating in a
# BsaI-digested PCR product.

setwd("/home/claude/goldengateR")
for (f in list.files("R", full.names = TRUE)) source(f)

cat("=== Workflow: PCR-into-gel-purified-backbone ===\n\n")

# === Step 1: Build a Type 3 entry vector ===
# Realistic-ish layout: CamR + ColE1 (placeholder sequences) + BsaI-flanked
# GFP dropout with CCAT/AATG overhangs (the canonical YTK Type 3 overhangs
# marking the start and end of a CDS).
camR <- paste(rep("ACGT", 200), collapse = "")    # 800 nt placeholder for CamR
ori  <- paste(rep("CGAT", 150), collapse = "")    # 600 nt placeholder for ColE1 ori
gfp_dropout <- paste(rep("AAGGCC", 120), collapse = "")  # 720 nt placeholder for sfGFP

# Vector layout: sites face inward toward the GFP dropout, recognition
# stays with the dropout, backbone is clean.
entry_vec_seq <- paste0(
  camR,                       # 1..800
  "CCAT", "A", "GAGACC",      # 801..811: rev BsaI site, OH=CCAT (becomes backbone.right_oh)
  gfp_dropout,                # 812..1531: dropout (will get filtered/discarded in one-pot,
                              #             but here we'll gel-select the OTHER fragment)
  "GGTCTC", "A", "AATG",      # 1532..1542: fwd BsaI site, OH=AATG (becomes backbone.left_oh)
  ori                         # 1543..2142
)
entry_vec <- new_part(entry_vec_seq, topology = "circular", name = "pYTK_T3_dest")
cat(sprintf("Entry vector: %d bp circular\n", nchar(entry_vec_seq)))

# === Step 2: Digest and gel-purify the LARGE backbone fragment ===
backbone <- digest_and_select(entry_vec, enzyme = "BsaI", select = "largest")
cat(sprintf("\nBackbone after gel purification:\n"))
cat(sprintf("  Topology:  %s\n", backbone$topology))
cat(sprintf("  Length:    %d bp\n", nchar(backbone$sequence)))
cat(sprintf("  left_oh:   %s\n", attr(backbone, "left_oh")))
cat(sprintf("  right_oh:  %s\n", attr(backbone, "right_oh")))
stopifnot(attr(backbone, "left_oh") == "AATG")
stopifnot(attr(backbone, "right_oh") == "CCAT")
cat("PASS\n")

# === Step 3: Build the PCR product ===
# Backbone(L=AATG, R=CCAT) means in the assembly cycle:
#   backbone.right(CCAT) -> insert.left(must be CCAT)
#   insert.right(must be AATG) -> backbone.left(AATG)
# So the released insert needs left_oh=CCAT, right_oh=AATG.
#
# Primer design:
# Forward primer:  5'-CGAGCG-GGTCTC-A-CCAT-[anneal]-3'    (gives insert.left_oh = CCAT)
# Reverse primer:  5'-CGAGCG-GGTCTC-A-CATT-[anneal]-3'    (CATT = revcomp of AATG)
# After PCR, top-strand:
#   CGAGCG-GGTCTC-A-CCAT-[gene]-AATG-T-GAGACC-CGCTCG
gene_with_internal_site <- paste0(
  "ATGAAACCGCTGGAGAAACTGGCGAAATATGGCAGCATTGAA",
  "CTGGGTAACGGCAAATTAGCGAAAGAATGGCTGTAACCGCATAGGTACGTAGCTAGCATG",
  "AAACCGCTGGAGAAACTGGCGAAATATGGCAGCATTGAACTGGGTAACGGCAAATTAGCGAAAGAATGGCTGTAA"
)
forward_tail <- paste0("CGAGCG", "GGTCTC", "A", "CCAT")
reverse_tail_top <- paste0("AATG", "T", "GAGACC", "CGCTCG")
pcr_product_seq <- paste0(forward_tail, gene_with_internal_site, reverse_tail_top)
cat(sprintf("\nPCR product: %d bp linear\n", nchar(pcr_product_seq)))

writeLines(c(">xylB_amp", pcr_product_seq), "/tmp/xylB_amp.fa")
amp <- read_part("/tmp/xylB_amp.fa", topology = "linear")

# === Step 4: Digest the PCR product, take the middle fragment ===
insert <- digest_and_select(amp, enzyme = "BsaI", select = "largest")
cat(sprintf("\nInsert after digest:\n"))
cat(sprintf("  Length:    %d bp\n", nchar(insert$sequence)))
cat(sprintf("  left_oh:   %s\n", attr(insert, "left_oh")))
cat(sprintf("  right_oh:  %s\n", attr(insert, "right_oh")))
stopifnot(attr(insert, "left_oh") == "CCAT")
stopifnot(attr(insert, "right_oh") == "AATG")
cat("PASS\n")

# === Step 5: Ligate ===
asm <- ligate(list(backbone, insert), enzyme = "BsaI",
              name = "pYTK_T3_xylB",
              output = "/tmp/pYTK_T3_xylB.gb")
cat(sprintf("\nLigation product: %d bp circular\n", nchar(asm$sequence)))

# Sanity checks
stopifnot(asm$topology == "circular")
stopifnot(grepl(gene_with_internal_site, asm$sequence))   # insert is in
stopifnot(!grepl(gfp_dropout, asm$sequence))               # dropout is out
stopifnot(grepl(camR, asm$sequence))                       # CamR retained
stopifnot(!grepl("GGTCTC", asm$sequence))                  # sites consumed
stopifnot(!grepl("GAGACC", asm$sequence))
cat("All sanity checks PASS\n")

# === Step 6: Round-trip the GenBank ===
back <- read_part("/tmp/pYTK_T3_xylB.gb")
stopifnot(back$sequence == asm$sequence)
cat("Round-trip PASS\n")

# === Step 7: Test selection by size band ===
# The entry vector digest gives a 1408 bp backbone and a 742 bp dropout.
bb_band <- digest_and_select(entry_vec, "BsaI", select = c(1300, 1500))
stopifnot(nchar(bb_band$sequence) == nchar(backbone$sequence))
cat("Size-band selection (1300-1500 bp): PASS\n")

# Out-of-range selection should error
err <- tryCatch(
  digest_and_select(entry_vec, "BsaI", select = c(5000, 10000)),
  error = function(e) e
)
stopifnot(inherits(err, "error"))
cat("Out-of-range selection correctly errors: PASS\n")

# Ambiguous size band (both fragments fall in range) should error
err <- tryCatch(
  digest_and_select(entry_vec, "BsaI", select = c(0, 10000)),
  error = function(e) e
)
stopifnot(inherits(err, "error"))
cat("Ambiguous size band correctly errors: PASS\n")

# === Step 8: make_fragment for synthetic / oligo input ===
# Backbone has L=AATG, R=CCAT, so oligo needs L=CCAT, R=AATG to chain.
oligo_frag <- make_fragment(left_oh = "CCAT",
                            sequence = "ATGAAATAA",
                            right_oh = "AATG",
                            name = "short_oligo")
stopifnot(attr(oligo_frag, "left_oh") == "CCAT")
stopifnot(nchar(oligo_frag$sequence) == 4 + 9 + 4)
asm2 <- ligate(list(backbone, oligo_frag), enzyme = "BsaI",
               name = "oligo_construct")
stopifnot(grepl("ATGAAATAA", asm2$sequence))
cat("Synthetic oligo via make_fragment(): PASS\n")

# === Step 9: Wrong-overhang fragment should fail to ligate ===
bad_oligo <- make_fragment(left_oh = "TTTT", sequence = "ATG",
                           right_oh = "AAAA", name = "wrong_oh")
err <- tryCatch(
  ligate(list(backbone, bad_oligo)),
  error = function(e) e
)
stopifnot(inherits(err, "error"))
cat("Mismatched overhangs correctly fail: PASS\n")

cat("\n========== ALL WORKFLOW TESTS PASSED ==========\n")
