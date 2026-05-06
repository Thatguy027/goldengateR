setwd("/home/claude/goldengateR")
for (f in list.files("R", full.names = TRUE)) source(f)

# Build a vector and insert with annotated features, write to GenBank,
# read back, assemble, and verify the features carried through.

ovhg1 <- "TGCC"; ovhg2 <- "ATGG"
left_arm  <- "AAAAAAAAAACCCCCCCCCC"
right_arm <- "GGGGGGGGGGGGGGGGGGGG"
dropout   <- "TTTTTTTTTTTTTTTTTTTT"
entry_arm <- "CACACACACACACACACACA"
insert_core <- "ATATATATATATATATATAT"

vec_seq <- paste0(left_arm, ovhg1, "A", "GAGACC", dropout, "GGTCTC", "A", ovhg2, right_arm)
vec_features <- data.frame(
  type = c("misc_feature", "CDS"),
  start = c(1L, 32L),
  end = c(20L, 51L),
  strand = c(1L, 1L),
  stringsAsFactors = FALSE
)
vec_features$qualifiers <- list(
  list(label = "left_arm_marker", note = "vector left arm"),
  list(label = "dropout_GFP", product = "drop-out marker")
)
vec_part <- new_part(vec_seq, features = vec_features, topology = "circular", name = "test_vec")

ins_seq <- paste0(entry_arm, "GGTCTC", "A", ovhg1, insert_core, ovhg2, "A", "GAGACC", entry_arm)
# Insert positions: entry_arm 1-20, GGTCTC 21-26, A 27, ovhg1 28-31, insert_core 32-51, ovhg2 52-55
ins_features <- data.frame(
  type = c("CDS"),
  start = c(32L),
  end = c(51L),
  strand = c(1L),
  stringsAsFactors = FALSE
)
ins_features$qualifiers <- list(list(label = "my_gene", product = "test product"))
ins_part <- new_part(ins_seq, features = ins_features, topology = "circular", name = "test_ins")

# Write both parts to GenBank
write_genbank(vec_part, "/tmp/test_vec.gb")
write_genbank(ins_part, "/tmp/test_ins.gb")
cat("--- /tmp/test_vec.gb ---\n")
cat(paste(readLines("/tmp/test_vec.gb"), collapse = "\n"), "\n")
cat("\n--- /tmp/test_ins.gb ---\n")
cat(paste(readLines("/tmp/test_ins.gb"), collapse = "\n"), "\n")

# Read back from disk and assemble
v2 <- read_part("/tmp/test_vec.gb")
i2 <- read_part("/tmp/test_ins.gb")

cat("\nVec features after round-trip:\n")
print(v2$features[, c("type","start","end","strand")])
cat("Ins features after round-trip:\n")
print(i2$features[, c("type","start","end","strand")])

asm <- golden_gate_assemble(list(v2, i2), enzyme = "BsaI",
                            name = "annotated_asm",
                            output = "/tmp/annotated_asm.gb")

cat("\nAssembly features:\n")
print(asm$features[, c("type","start","end","strand")])
cat("Assembly qualifiers:\n")
for (i in seq_len(nrow(asm$features))) {
  q <- asm$features$qualifiers[[i]]
  cat(sprintf("  Feature %d (%s, %d..%d): %s\n", i, asm$features$type[i],
              asm$features$start[i], asm$features$end[i],
              paste(names(q), unlist(q), sep="=", collapse=", ")))
}

# Expected: the my_gene CDS from the insert should be present, dropout_GFP
# should NOT be present (it was in the dropout, which got filtered out).
qualifiers_flat <- unlist(lapply(asm$features$qualifiers, function(q) {
  if (is.null(q)) return(character())
  vapply(seq_along(q), function(i) paste(names(q)[i], q[[i]], sep="="), character(1))
}))
cat("\nAll qualifier values:\n", paste(qualifiers_flat, collapse="\n  "), "\n")

stopifnot(any(grepl("my_gene", qualifiers_flat)))
stopifnot(any(grepl("left_arm_marker", qualifiers_flat)))
stopifnot(!any(grepl("dropout_GFP", qualifiers_flat)))
cat("\nFeature carry-through: PASS\n")

# Look at the assembly GenBank output
cat("\n--- /tmp/annotated_asm.gb ---\n")
cat(paste(readLines("/tmp/annotated_asm.gb"), collapse = "\n"), "\n")
