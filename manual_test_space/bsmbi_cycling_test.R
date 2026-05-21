# Test: BsmBI cycling reaction on xylb_cassette + pMYT039
# Expects 4 possible circular products:
#   - xylb re-circularized (original)
#   - pMYT039 re-circularized (original)
#   - two shuffled/recombined products

pkg_root <- rprojroot::find_package_root_file()
devtools::load_all(pkg_root, quiet = TRUE)

plasmid_dir <- file.path(pkg_root, "inst/extdata/plasmids")
out_dir     <- file.path(pkg_root, "manual_test_space/cycling_products")
dir.create(out_dir, showWarnings = FALSE)

xylb   <- read_part(file.path(plasmid_dir, "xylb_cassette.gb"))
pmyt39 <- read_part(file.path(plasmid_dir, "pMYT039.gb"))

enzyme <- "BsmBI"

cat("=== Digesting with", enzyme, "===\n\n")
frags_xylb <- goldengateR:::digest_part_no_filter(xylb,   enzyme)
frags_pmyt <- goldengateR:::digest_part_no_filter(pmyt39, enzyme)

print_frags <- function(name, frags) {
  cat(sprintf("%s (%d fragment(s)):\n", name, length(frags)))
  for (i in seq_along(frags)) {
    f <- frags[[i]]
    cat(sprintf("  [%d] left_oh=%-4s  right_oh=%-4s  core=%d bp\n",
                i, f$left_oh, f$right_oh, nchar(f$sequence)))
  }
  cat("\n")
}
print_frags("xylb_cassette", frags_xylb)
print_frags("pMYT039",       frags_pmyt)

all_frags <- c(frags_xylb, frags_pmyt)
n <- length(all_frags)
labels <- c(
  paste0("xylb[", seq_along(frags_xylb), "]"),
  paste0("pmyt39[", seq_along(frags_pmyt), "]")
)

cat("=== Trying all 2-fragment circular assemblies ===\n\n")
products <- list()

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    pair   <- all_frags[c(i, j)]
    result <- tryCatch(
      goldengateR:::find_circular_assembly(pair),
      error = function(e) NULL
    )
    if (is.null(result)) next

    stitched <- goldengateR:::stitch_assembly(result$path, result$oriented)
    pname <- paste0(labels[i], "_x_", labels[j])
    product <- new_part(
      sequence = stitched$sequence,
      features = stitched$features,
      topology = "circular",
      name     = pname
    )
    products[[length(products) + 1]] <- product

    out_file <- file.path(out_dir, paste0(pname, ".gb"))
    write_genbank(product, out_file)
    cat(sprintf("VALID: %s + %s  =>  %d bp  [%s]\n",
                labels[i], labels[j], nchar(product$sequence), out_file))
  }
}

cat(sprintf("\n=== %d circular product(s) assembled ===\n", length(products)))
if (length(products) == 0) {
  cat("No valid assemblies found — overhangs may not be compatible.\n")
}
