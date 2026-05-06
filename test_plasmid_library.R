# Test the bundled plasmid library and user-dir override.

setwd("/home/claude/goldengateR")
for (f in list.files("R", full.names = TRUE)) source(f)

cat("=== Test 1: list_plasmids() finds bundled plasmids ===\n")
# In dev mode (sourcing R/), the search-path fallback should locate
# inst/extdata/plasmids relative to the working dir.
all <- list_plasmids()
print(all)
stopifnot(nrow(all) >= 2)
stopifnot("pYTK_T3_dest_demo" %in% all$name)
stopifnot("pYTK_T3_xylB_demo" %in% all$name)
stopifnot(all(all$source == "bundled"))
cat("PASS\n\n")

cat("=== Test 2: load_plasmid() by name ===\n")
dest <- load_plasmid("pYTK_T3_dest_demo")
stopifnot(inherits(dest, "ggr_part"))
stopifnot(dest$topology == "circular")
stopifnot(nchar(dest$sequence) == 2324)
cat(sprintf("Loaded %s: %d bp circular\n", dest$name, nchar(dest$sequence)))
cat("PASS\n\n")

cat("=== Test 3: case-insensitive lookup ===\n")
dest2 <- load_plasmid("PYTK_T3_DEST_DEMO")
stopifnot(dest2$sequence == dest$sequence)
cat("PASS\n\n")

cat("=== Test 4: helpful error on missing plasmid ===\n")
err <- tryCatch(load_plasmid("does_not_exist"), error = function(e) e)
stopifnot(inherits(err, "error"))
cat("Error message:", conditionMessage(err), "\n")
stopifnot(grepl("does_not_exist", conditionMessage(err)))
stopifnot(grepl("Available", conditionMessage(err)))
cat("PASS\n\n")

cat("=== Test 5: pattern filter ===\n")
matches <- list_plasmids(pattern = "xylB")
stopifnot(nrow(matches) == 1)
stopifnot(matches$name == "pYTK_T3_xylB_demo")
cat("PASS\n\n")

cat("=== Test 6: user-configured plasmid directory takes precedence ===\n")
# Create a user dir with a same-named file and verify it overrides bundled
user_dir <- tempfile("user_plasmids_")
dir.create(user_dir)
# Write a different plasmid with the SAME name as a bundled one
fake <- new_part("ATCG", topology = "linear", name = "pYTK_T3_dest_demo")
write_genbank(fake, file.path(user_dir, "pYTK_T3_dest_demo.gb"))
set_plasmid_dir(user_dir)

available <- list_plasmids()
print(available)
# The bundled "pYTK_T3_dest_demo" should be shadowed by the user version
override <- available[available$name == "pYTK_T3_dest_demo", ]
stopifnot(nrow(override) == 1)
stopifnot(override$source == "user")

# Loading by name should give the user version (which has only "ATCG" sequence)
loaded_user <- load_plasmid("pYTK_T3_dest_demo")
stopifnot(loaded_user$sequence == "ATCG")
cat("User dir override: PASS\n")

# The bundled-only plasmid should still be accessible
bundled_only <- load_plasmid("pYTK_T3_xylB_demo")
stopifnot(nchar(bundled_only$sequence) == 1500)
cat("Bundled fallback when not shadowed: PASS\n")

# Clear and verify back to bundled-only
set_plasmid_dir(NULL)
back_to_bundled <- list_plasmids()
stopifnot(all(back_to_bundled$source == "bundled"))
stopifnot(nrow(back_to_bundled) == 2)
cat("Clearing user dir: PASS\n\n")

cat("=== Test 7: end-to-end - use load_plasmid() in an assembly ===\n")
dest <- load_plasmid("pYTK_T3_dest_demo")
backbone <- digest_and_select(dest, "BsaI", select = "largest")

amp_seq <- paste0("CGAGCG", "GGTCTC", "A", "AATG",
                  "ATGAAATTTAACGGCTGCTAA",
                  "CCAT", "T", "GAGACC", "CGCTCG")
amp <- new_part(amp_seq, topology = "linear", name = "test_gene_amp")
insert <- digest_and_select(amp, "BsaI", select = "largest")

asm <- suppressWarnings(ligate(list(backbone, insert), enzyme = "BsaI",
                                name = "test_construct"))
stopifnot(asm$topology == "circular")
stopifnot(grepl("ATGAAATTTAACGGCTGCTAA", asm$sequence))
cat(sprintf("Assembled %d bp construct via load_plasmid: PASS\n",
            nchar(asm$sequence)))

cat("\n========== ALL PLASMID LIBRARY TESTS PASSED ==========\n")
