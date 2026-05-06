# Test script - run with: Rscript test_run.R
# Loads all R/ files manually (no formal install needed) and runs sanity checks
# against hand-computed expected outcomes.

setwd("/home/claude/goldengateR")
for (f in list.files("R", full.names = TRUE)) source(f)

cat("=== Test 1: revcomp ===\n")
stopifnot(revcomp("GGTCTC") == "GAGACC")
stopifnot(revcomp("AAAA") == "TTTT")
stopifnot(revcomp("ACGT") == "ACGT")
cat("PASS\n\n")

cat("=== Test 2: site finding on a designed substrate ===\n")
# Construct a circular plasmid with two BsaI sites in inverted orientation
# (canonical YTK Type 0 entry vector layout): the two sites flank a dropout
# such that when cut, they release the dropout and leave the backbone with
# overhangs that will accept an insert.
#
# Design:
#   ...AAAA GGTCTC A CCCC [DROPOUT] AAAA T GAGACC TTTT...
# Forward BsaI site GGTCTC at pos X cuts top after X+6 (= A spacer end), with
# overhang = the 4 bases after it ("CCCC").
# Reverse BsaI site GAGACC at pos Y cuts top 5 nt UPSTREAM of Y (i.e. removes
# 4 nt overhang on the downstream side). Overhang here is the 4 nt
# immediately preceding GAGACC on the top strand (TTTT? no — the 4 nt
# starting at Y - 5).
#
# Build a clean test:
backbone <- paste0(
  "AAAAAAAAAA",            # 1-10  (left arm)
  "TTTT",                  # 11-14 left overhang (will appear at junction)
  "A",                     # 15    spacer
  "GAGACC",                # 16-21 reverse BsaI (cuts upstream)
  "GGGGGGGGGGGGGGGG",      # 22-37 backbone middle
  "GGTCTC",                # 38-43 forward BsaI
  "A",                     # 44    spacer
  "CCCC",                  # 45-48 right overhang
  "AAAAAAAAAA"             # 49-58 (right arm)
)
# After BsaI digestion, the backbone fragment has:
#   - left overhang TTTT (from the reverse site)
#   - right overhang CCCC (from the forward site)
#   - the recognition+spacer is excised on each end
#
# Forward site at 38: top cut after 38+5+1=44, bottom cut after 38+5+5=48.
# Overhang = positions 45..48 = "CCCC". Backbone right end keeps 1..44.
# Wait — that's wrong, the forward cut RELEASES the backbone INTO the right
# arm. Let me rethink for circular topology.
#
# In a circular plasmid with the layout above, BsaI cuts produce two
# fragments: the "dropout" (between the two cuts going one way) and the
# "backbone" (going the other way). The fragment that contains positions
# 22-37 (the GGGG middle) is the dropout because the two enzyme sites face
# OUTWARD from it. To make the GGGG the backbone instead, the sites must
# face INWARD toward the dropout.

# Redesign with inward-facing sites flanking a dropout:
#   left arm | GGTCTC N OVHG1 [DROPOUT] OVHG2 N GAGACC | right arm
# Both sites cut INTO the dropout, releasing it. Backbone retains the
# regions outside the cut points and joins them via the circular topology.
# After digestion, the backbone fragment (which includes the left+right arms
# joined through the origin) has:
#   - left overhang (on the side facing where the dropout used to be on the
#     right of the right arm) = OVHG2
#   - right overhang (on the side facing where the dropout used to be on the
#     left of the left arm) = OVHG1
# Wait, I'm getting confused by geometry. Let me just build it and inspect.

cat("Building synthetic plasmid 1 (vector with dropout)...\n")
left_arm  <- "AAAAAAAAAACCCCCCCCCC"          # 20 nt
ovhg1     <- "TGCC"                          # left overhang (top-strand)
ovhg2     <- "ATGG"                          # right overhang
right_arm <- "GGGGGGGGGGGGGGGGGGGG"          # 20 nt
dropout   <- "TTTTTTTTTTTTTTTTTTTT"          # 20 nt
# Vector layout (circular):
#   left_arm + GAGACC + N + ovhg1 + dropout + ovhg2 + N + GGTCTC + right_arm
# Reading left to right in the linear representation:
# - Reverse site GAGACC just after left_arm: cuts UPSTREAM on top strand, so
#   it cuts inside left_arm? No - reverse-strand BsaI at position p cuts top
#   strand at position p - 5 (between p-5 and p-4). So if GAGACC is at pos
#   21..26, top cut is after pos 16. The 4 nt overhang on the downstream
#   fragment (containing the dropout) starts at pos 17.
# This means I have to put OVHG1 BEFORE the GAGACC, not after it.
#
# Cleaner design:
#   left_arm + ovhg1 + N + GAGACC + dropout_core + GGTCTC + N + ovhg2 + right_arm
# where the dropout_core has no enzyme sites and is what gets released.
plasmid1_seq <- paste0(
  left_arm,    # 1-20
  ovhg1,       # 21-24  becomes left overhang of dropout, right overhang of backbone (in top-strand coords at the junction = ovhg1)
  "A",         # 25     spacer N (consumed by enzyme cut, ends up in discarded fragment)
  "GAGACC",    # 26-31  reverse BsaI site
  dropout,     # 32-51  20 nt dropout core
  "GGTCTC",    # 52-57  forward BsaI site
  "A",         # 58     spacer N
  ovhg2,       # 59-62  becomes right overhang of dropout, left overhang of backbone (top-strand at junction = ovhg2)
  right_arm    # 63-82
)
nchar(plasmid1_seq)  # 82

# Verify cut positions on this:
sites <- find_sites(plasmid1_seq, "BsaI", circular = TRUE)
print(sites)
# Expect:
# Forward site GGTCTC at pos 52: top cut after 52+6+1-1 = 58, bottom cut after 52+6+5-1 = 62
#   overhang at positions 59..62 = ovhg2 = "ATGG"  ✓
# Reverse site GAGACC at pos 26: top cut after 26 - 5 - 1 = 20, bottom cut after 26 - 1 - 1 = 24
#   overhang at positions 21..24 = ovhg1 = "TGCC"  ✓
stopifnot(nrow(sites) == 2)
oh1 <- substring(plasmid1_seq, sites$oh_start[1], sites$oh_end[1])
oh2 <- substring(plasmid1_seq, sites$oh_start[2], sites$oh_end[2])
cat("Found overhangs:", oh1, "and", oh2, "\n")
stopifnot(setequal(c(oh1, oh2), c(ovhg1, ovhg2)))
cat("PASS\n\n")

cat("=== Test 3: digest the vector ===\n")
vector_part <- new_part(plasmid1_seq, topology = "circular", name = "vec1")
frags <- digest_part(vector_part, "BsaI")
cat("Got", length(frags), "fragments after digestion\n")
for (i in seq_along(frags)) {
  cat(sprintf("  Frag %d: left_oh=%s, right_oh=%s, core_len=%d, core_preview=%s\n",
              i, frags[[i]]$left_oh, frags[[i]]$right_oh,
              nchar(frags[[i]]$sequence),
              substr(frags[[i]]$sequence, 1, 20)))
}
# Expect 1 usable fragment (the backbone). The dropout fragment carries the
# recognition sites (they face inward toward the dropout), so it is
# correctly filtered out as "would be re-cut".
stopifnot(length(frags) == 1)
stopifnot(frags[[1]]$left_oh == ovhg2)
stopifnot(frags[[1]]$right_oh == ovhg1)
cat("PASS\n\n")

cat("=== Test 4: build a matching insert and assemble ===\n")
# Build an insert plasmid (Level 0 entry vector containing the new insert
# flanked by BsaI sites that release it with overhangs ovhg1 and ovhg2,
# matching the vector's overhangs so the insert ligates into the backbone).
#
# Layout for entry vector with insert that goes BACKBONE -> insert -> BACKBONE:
# The released insert fragment should have:
#   left overhang = ovhg1   (matches the backbone's right overhang)
#   right overhang = ovhg2  (matches the backbone's left overhang)
# Wait — in a circular ligation, fragment A's right_oh must equal fragment
# B's left_oh. Backbone's right_oh comes from the forward-site cut and is
# ovhg2. So insert's left_oh must be ovhg2. Backbone's left_oh = ovhg1, so
# insert's right_oh must be ovhg1. Therefore insert layout:
#   entry_arm + GGTCTC + N + ovhg2 + INSERT_CORE + ovhg1 + N + GAGACC + entry_arm

entry_arm <- "CACACACACACACACACACA"  # 20 nt
insert_core <- "ATATATATATATATATATAT"  # 20 nt
# Build an insert plasmid (Level 0 entry vector). For the insert to ligate
# into the vector backbone, its overhangs must be the COMPLEMENT of the
# vector's: backbone left_oh = ATGG (from forward cut) so insert right_oh
# must = ATGG; backbone right_oh = TGCC (from reverse cut) so insert left_oh
# must = TGCC.
#
# To get insert left_oh = TGCC from a cut, we need a REVERSE site upstream
# of the insert (whose cut leaves TGCC as the overhang on the downstream
# fragment). To get insert right_oh = ATGG from a cut, we need a FORWARD
# site downstream of the insert. So the insert plasmid layout is the
# OPPOSITE orientation of the vector:
#   entry_arm + GAGACC + N + TGCC + INSERT_CORE + ATGG + N + GGTCTC + entry_arm
# Wait - reverse-site cut overhang is the 4 nt UPSTREAM of GAGACC, so
# we put TGCC immediately before GAGACC. And forward-site cut overhang is
# the 4 nt DOWNSTREAM of GGTCTC+spacer, but here we want it on the LEFT of
# GGTCTC (because we want it inside the insert, before the site). So the
# forward site must be flipped... no, "forward" just means the site reads
# GGTCTC on top strand. We need the overhang to come from the insert side.
#
# Actually the cleanest way: the insert plasmid layout has the recognition
# sites pointing OUTWARD (away from the insert, into the discarded entry
# arms) — opposite to the vector where they point INWARD toward the dropout.
#
# Outward-facing means: forward site (GGTCTC) is to the LEFT of the insert
# (cuts to its right, into the insert); reverse site (GAGACC) is to the
# RIGHT of the insert (cuts to its left, into the insert). Wait, that's
# inward-facing for the insert, same as the vector... let me just compute.
#
# Forward GGTCTC + N + XXXX + ... : XXXX is the overhang on what's downstream
#   of the cut. If insert is downstream, insert's LEFT overhang = XXXX.
# Reverse: ... + YYYY + N + GAGACC : YYYY is overhang on what's UPSTREAM of
#   the cut (= upstream of GAGACC). If insert is upstream, insert's RIGHT
#   overhang = YYYY.
# So plasmid layout: entry_arm + GGTCTC + N + LEFT_OH + INSERT + RIGHT_OH + N + GAGACC + entry_arm
# is exactly the VECTOR layout but with the insert in place of the dropout.
# This produces an insert with left_oh = LEFT_OH and right_oh = RIGHT_OH.
#
# Setting LEFT_OH = TGCC and RIGHT_OH = ATGG gives an insert whose overhangs
# are swapped relative to my first attempt.
plasmid2_seq <- paste0(
  entry_arm,
  "GGTCTC",     # forward BsaI
  "A",
  ovhg1,        # = TGCC: this becomes left overhang of insert
  insert_core,
  ovhg2,        # = ATGG: this becomes right overhang of insert
  "A",
  "GAGACC",     # reverse BsaI
  entry_arm
)
nchar(plasmid2_seq)

insert_part <- new_part(plasmid2_seq, topology = "circular", name = "ent1")
ifrags <- digest_part(insert_part, "BsaI")
cat("Insert digest: got", length(ifrags), "fragments\n")
for (i in seq_along(ifrags)) {
  cat(sprintf("  Frag %d: left_oh=%s, right_oh=%s, core_len=%d\n",
              i, ifrags[[i]]$left_oh, ifrags[[i]]$right_oh,
              nchar(ifrags[[i]]$sequence)))
}

# Assembly
asm <- golden_gate_assemble(list(vector_part, insert_part),
                             enzyme = "BsaI",
                             name = "test_asm")
cat("\nAssembled product:\n")
print(asm)
cat("Sequence length:", nchar(asm$sequence), "bp\n")
cat("Sequence:\n", asm$sequence, "\n", sep = "")

# Sanity checks on the assembled product:
#  - Should be circular
#  - Should not contain GGTCTC or GAGACC (Type IIS sites consumed)
#  - Should contain the insert_core sequence
#  - Should contain left_arm and right_arm and entry_arm
stopifnot(asm$topology == "circular")
stopifnot(!grepl("GGTCTC", asm$sequence))
stopifnot(!grepl("GAGACC", asm$sequence))
stopifnot(grepl(insert_core, asm$sequence))
# The dropout sequence should NOT be in the product
stopifnot(!grepl(dropout, asm$sequence))
# Both arms of the vector backbone should be present (possibly across the origin)
double_seq <- paste0(asm$sequence, asm$sequence)
stopifnot(grepl(left_arm, double_seq))
stopifnot(grepl(right_arm, double_seq))
# entry_arm should NOT be in the product - it's part of the discarded entry vector fragment
stopifnot(!grepl(entry_arm, double_seq))
cat("\nPASS - all sanity checks on assembly product\n\n")

cat("=== Test 5: write and re-read GenBank ===\n")
write_genbank(asm, "/tmp/test_asm.gb")
cat("Wrote /tmp/test_asm.gb\n")
re_read <- read_part("/tmp/test_asm.gb")
cat("Re-read OK\n")
print(re_read)
stopifnot(re_read$sequence == asm$sequence)
stopifnot(re_read$topology == "circular")
cat("Round-trip identical: PASS\n\n")

cat("=== Test 6: ambiguous assembly should error ===\n")
# Build a third plasmid with the SAME overhangs as the insert -> two
# possible assemblies (vector + insert1, or vector + insert2)
plasmid3_seq <- paste0(
  entry_arm,
  "GGTCTC", "A", ovhg1,         # same overhangs as insert
  "GCGCGCGCGCGCGCGCGCGC",       # different core
  ovhg2, "A", "GAGACC",
  entry_arm
)
insert2_part <- new_part(plasmid3_seq, topology = "circular", name = "ent2")
err <- tryCatch(
  golden_gate_assemble(list(vector_part, insert_part, insert2_part),
                       enzyme = "BsaI"),
  error = function(e) e
)
stopifnot(inherits(err, "error"))
cat("Got expected error:", conditionMessage(err), "\n")
cat("PASS\n\n")

cat("=== Test 7: impossible assembly should error ===\n")
# Vector with overhangs that don't match the insert
mismatched_seq <- paste0(
  entry_arm,
  "GGTCTC", "A", "AAAA",  # ovhg = AAAA, not ovhg1 or ovhg2
  "GCGCGCGCGCGCGCGCGCGC",
  "TTTT", "A", "GAGACC",
  entry_arm
)
mismatched_part <- new_part(mismatched_seq, topology = "circular", name = "mm")
err <- tryCatch(
  golden_gate_assemble(list(vector_part, mismatched_part), enzyme = "BsaI"),
  error = function(e) e
)
stopifnot(inherits(err, "error"))
cat("Got expected error:", conditionMessage(err), "\n")
cat("PASS\n\n")

cat("=== Test 8: 3-part assembly (more realistic YTK-style) ===\n")
# Three-part circular assembly: vector backbone + part A + part B
# Cycle order: backbone(L=oh_a, R=oh_c) -> A(L=oh_c, R=oh_b) -> B(L=oh_b, R=oh_a) -> backbone
oh_a <- "AATG"
oh_b <- "GCTT"
oh_c <- "CGCT"
# Vector layout: sites face INWARD toward the dropout so recognition sites
# are retained in the (filtered-out) dropout fragment and the backbone is clean.
#   [left_arm][OH_BACKBONE_L][N][GAGACC][dropout][GGTCTC][N][OH_BACKBONE_R][right_arm]
# Reverse site upstream of dropout: its cut leaves OH on the upstream fragment's
#   right side (= backbone right_oh, via origin wrap).
# Forward site downstream of dropout: its cut leaves OH on the downstream
#   fragment's left side (= backbone left_oh, via origin wrap).
# Wait, let me re-derive once more from scratch using the plasmid1 result
# which I know works:
#   plasmid1: [arm1][TGCC][A][GAGACC][dropout][GGTCTC][A][ATGG][arm2]
#   digest gave backbone(L=ATGG, R=TGCC). So:
#     forward-site OH (ATGG, on the right in linear) -> backbone left_oh
#     reverse-site OH (TGCC, on the left in linear)  -> backbone right_oh
#   This is consistent with: backbone wraps origin, so its "left_oh" is
#   actually at the right end of arm2-region in linear coords (just past
#   GGTCTC+spacer), and its "right_oh" is at the left end of arm1-region
#   in linear coords (just before GAGACC).
# So to set backbone(L=oh_a, R=oh_c), I want:
#   forward-site OH (downstream of GGTCTC+spacer) = oh_a
#   reverse-site OH (upstream of GAGACC) = oh_c
# i.e.:
vec3_seq <- paste0(
  "AAAAAAAAAACCCCCCCCCC",                  # left_arm (arm1)
  oh_c, "A", "GAGACC",                     # reverse site, OH=oh_c on its left
  "GGGGGGGGGGGGGGGGGGGG",                  # dropout
  "GGTCTC", "A", oh_a,                     # forward site, OH=oh_a on its right
  "TTTTTTTTTTTTTTTTTTTT"                   # right_arm (arm2)
)
vec3 <- new_part(vec3_seq, topology = "circular", name = "vec3")
fr_v <- digest_part(vec3, "BsaI")
cat("vec3 digest gave", length(fr_v), "fragments:\n")
for (f in fr_v) cat("  L=", f$left_oh, " R=", f$right_oh, " core_len=", nchar(f$sequence), "\n", sep="")
# Use whatever the digest actually says to design A and B
backbone_L <- fr_v[[1]]$left_oh
backbone_R <- fr_v[[1]]$right_oh
cat("Backbone overhangs: L=", backbone_L, " R=", backbone_R, "\n", sep="")

# In the cycle backbone.R -> A.L, A.R -> B.L, B.R -> backbone.L
# So A.L = backbone_R, B.R = backbone_L. Pick A.R = B.L = oh_b (free choice).
A_L <- backbone_R
B_R <- backbone_L
A_R <- oh_b
B_L <- oh_b

# Part A plasmid (insert layout: fwd-N-OH_FWD-insert-OH_REV-N-rev gives
# released fragment with left_oh=OH_FWD, right_oh=OH_REV)
partA_seq <- paste0(
  "CACACACACACACACACACA",
  "GGTCTC", "A", A_L,
  "ATATATATATATATATATAT",
  A_R, "A", "GAGACC",
  "CACACACACACACACACACA"
)
partB_seq <- paste0(
  "TGTGTGTGTGTGTGTGTGTG",
  "GGTCTC", "A", B_L,
  "CCAACCAACCAACCAACCAA",
  B_R, "A", "GAGACC",
  "TGTGTGTGTGTGTGTGTGTG"
)

vec3 <- new_part(vec3_seq, topology = "circular", name = "vec3")
pA <- new_part(partA_seq, topology = "circular", name = "partA")
pB <- new_part(partB_seq, topology = "circular", name = "partB")

asm3 <- golden_gate_assemble(list(vec3, pA, pB), enzyme = "BsaI", name = "asm3")
print(asm3)
cat("3-part assembly length:", nchar(asm3$sequence), "bp\n")
stopifnot(grepl("ATATATATATATATATATAT", asm3$sequence))
stopifnot(grepl("CCAACCAACCAACCAACCAA", asm3$sequence))
stopifnot(!grepl("GGTCTC", asm3$sequence))
stopifnot(!grepl("GAGACC", asm3$sequence))
# Dropout should not be present
stopifnot(!grepl("GGGGGGGGGGGGGGGGGGGG", asm3$sequence))
cat("PASS\n\n")

cat("=== Test 9: BsmBI works the same way ===\n")
# Replace BsaI sites with BsmBI sites in plasmid1
plasmid1_bsmbi <- gsub("GGTCTC", "CGTCTC", plasmid1_seq)
plasmid1_bsmbi <- gsub("GAGACC", "GAGACG", plasmid1_bsmbi)  # revcomp(CGTCTC) = GAGACG
plasmid2_bsmbi <- gsub("GGTCTC", "CGTCTC", plasmid2_seq)
plasmid2_bsmbi <- gsub("GAGACC", "GAGACG", plasmid2_bsmbi)
v_b <- new_part(plasmid1_bsmbi, topology = "circular", name = "vec_b")
i_b <- new_part(plasmid2_bsmbi, topology = "circular", name = "ins_b")
asm_b <- golden_gate_assemble(list(v_b, i_b), enzyme = "BsmBI", name = "bsmbi_asm")
stopifnot(grepl(insert_core, asm_b$sequence))
stopifnot(!grepl("CGTCTC", asm_b$sequence))
cat("PASS\n\n")

cat("\n========== ALL TESTS PASSED ==========\n")
