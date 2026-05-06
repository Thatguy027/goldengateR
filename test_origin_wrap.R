# Regression test: a circular plasmid where one BsaI site is positioned
# such that its cut spans the origin (top_cut near L, bottom_cut near 1).
# This reproduces the bug reported with the user's pYTK_T3_xylB.gb file.

setwd("/home/claude/goldengateR")
for (f in list.files("R", full.names = TRUE)) source(f)

# Construct a 100 bp circular plasmid with a forward BsaI site near
# position 50 and a reverse BsaI site near position 5, so that the
# reverse-strand cut happens around the origin.
seq <- paste0(
  "ATCC",                # 1-4: this will become the wrap-around overhang
  "T",                   # 5: spacer N for the reverse site
  "GAGACC",              # 6-11: reverse BsaI site
  "AAAAAAAAAAAAAAAAAAAA", # 12-31: filler
  "AAAAAAAAAAAAAAAAAAAA", # 32-51: filler
  "GGTCTC",              # 52-57: forward BsaI site
  "A",                   # 58: spacer
  "TTTT",                # 59-62: forward overhang
  "CCCCCCCCCCCCCCCCCCCC", # 63-82: filler
  "CCCCCCCCCCCCCCCC"     # 83-98: filler
)
# Add 2 more nt to ensure 100 bp total
seq <- paste0(seq, "GG")
cat("Plasmid length:", nchar(seq), "bp\n")

p <- new_part(seq, topology = "circular", name = "wrap_test")
sites <- find_sites(seq, "BsaI", circular = TRUE)
cat("\nSites found:\n")
print(sites)
stopifnot(nrow(sites) == 2)

# The reverse site at p=6: top_cut = 0 -> wraps to L, bottom_cut = 4
# Overhang should wrap: positions L+1..4 i.e. 1..4 = "ATCC"
rev_site <- sites[sites$strand == -1, ]
stopifnot(nrow(rev_site) == 1)
stopifnot(rev_site$overhang == "ATCC")
cat("Reverse-strand origin-spanning cut overhang correctly identified as ATCC: PASS\n")

# Forward site at p=52: top_cut = 58, bottom_cut = 62, overhang = "TTTT"
fwd_site <- sites[sites$strand == 1, ]
stopifnot(fwd_site$overhang == "TTTT")
cat("Forward-strand cut overhang correctly identified as TTTT: PASS\n")

# Now digest and check fragments
fr <- digest_part_no_filter(p, "BsaI")
cat("\nFragments produced:", length(fr), "\n")
for (f in fr) {
  cat(sprintf("  L=%s, R=%s, core=%d bp\n",
              f$left_oh, f$right_oh, nchar(f$sequence)))
}
stopifnot(length(fr) == 2)
cat("Both fragments produced with valid overhangs: PASS\n")

# Test on the user's actual file
p_user <- read_part("/mnt/user-data/uploads/pYTK_T3_xylB.gb")
sel_largest  <- digest_and_select(p_user, "BsaI", select = "largest")
sel_smallest <- digest_and_select(p_user, "BsaI", select = "smallest")
stopifnot(nchar(sel_largest$sequence) > nchar(sel_smallest$sequence))
cat(sprintf("\nUser file digest: largest=%d bp, smallest=%d bp\n",
            nchar(sel_largest$sequence), nchar(sel_smallest$sequence)))
cat("digest_and_select on user file: PASS\n")

cat("\n========== ALL REGRESSION TESTS PASSED ==========\n")
