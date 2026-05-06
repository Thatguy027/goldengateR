# goldengateR

In silico Golden Gate assembly with BsaI and BsmBI in R, with annotated GenBank output.

## What it does

Reads plasmids and PCR products from GenBank or FASTA, simulates digestion with BsaI or BsmBI on both strands, ligates the resulting fragments via 4 nt sticky-end matching, and returns the unique circular product as an annotated GenBank file. Errors if the assembly is impossible or ambiguous.

Two assembly modes are supported. **One-pot** (`golden_gate_assemble`) digests and ligates everything in a single call — the right tool for standard YTK-style workflows. **Sequential** (`digest_and_select` + `ligate`) models the case where you cut a vector, gel-purify the backbone band, and ligate it with a separately-prepared PCR product.

A small library of bundled example plasmids ships in `inst/extdata/plasmids/` and is accessible by name via `load_plasmid()`. You can add your own files to that directory or point at a personal collection elsewhere.

## Install

```r
devtools::install_github("Thatguy027/goldengateR")
library(goldengateR)
```

No Bioconductor dependencies. After installation, `?goldengateR` shows the package overview and `?function_name` works for each exported function.

## Function reference

| Function | Purpose |
|---|---|
| `read_part(file, topology)` | Load a GenBank or FASTA file by path |
| `write_genbank(part, file)` | Write a `ggr_part` to GenBank |
| `load_plasmid(name)` | Load a plasmid from the library by name |
| `list_plasmids(pattern)` | See what plasmids are available |
| `set_plasmid_dir(dir)` | Add a supplemental plasmid directory |
| `plasmid_search_path()` | Inspect which directories are searched |
| `golden_gate_assemble(parts, enzyme, name, output)` | One-pot digest + ligate of circular plasmids |
| `digest_and_select(part, enzyme, select)` | Digest, then return one fragment by gel-band size |
| `make_fragment(left_oh, sequence, right_oh, name)` | Construct a fragment from sequence + overhangs |
| `ligate(fragments, enzyme, name, output)` | Ligate pre-prepared fragments into a circular product |
| `assemble_interactive(enzyme, output_dir)` | Prompt-driven single assembly at the console |
| `assemble_from_table(file, output_dir, enzyme)` | Batch assembly from a TSV/CSV manifest |

Full documentation for each function is available via `?function_name` after `library(goldengateR)`.

## Quick examples

### One-pot Golden Gate

```r
vec  <- load_plasmid("pYTK090")
p1   <- load_plasmid("pYTK008")
p234 <- load_plasmid("pYTK047")

asm <- golden_gate_assemble(
  parts   = list(vec, p1, p234),
  enzyme  = "BsaI",
  name    = "pYTK096_rebuilt",
  output  = "pYTK096_rebuilt.gb"
)
```

### Sequential digest + gel purification + ligation

```r
# 1. Digest the entry vector and gel-purify the large backbone band
entry_vec <- load_plasmid("pYTK_T3_dest_demo")
backbone  <- digest_and_select(entry_vec, "BsaI", select = "largest")

# 2. PCR-amplify your gene with BsaI primer tails, digest, take the insert
amp    <- read_part("xylB_amplicon.fa", topology = "linear")
insert <- digest_and_select(amp, "BsaI", select = "largest")

# 3. Ligate
asm <- ligate(list(backbone, insert),
              enzyme = "BsaI",
              name   = "pYTK_T3_xylB",
              output = "pYTK_T3_xylB.gb")
```

The key difference from `golden_gate_assemble`: gel-purified fragments that contain internal recognition sites are not filtered out, because in the wet-lab workflow the gel separates the band from the enzyme.

### Synthetic / annealed-oligo fragments

```r
oligo <- make_fragment(left_oh  = "CCAT",
                       sequence = "ATGAAAGCT",
                       right_oh = "AATG",
                       name     = "short_cds")
asm <- ligate(list(backbone, oligo), name = "oligo_construct")
```

## Plasmid library

```r
list_plasmids()                         # everything visible
list_plasmids(pattern = "pYTK0")        # filter by regex
load_plasmid("pYTK_T3_dest_demo")       # load by name (case-insensitive)

set_plasmid_dir("~/lab_plasmids/")      # add a supplemental directory
plasmid_search_path()                   # see search order
```

The bundled `inst/extdata/plasmids/` directory (pYTK001–096 plus demo plasmids) is always searched first. A directory set via `set_plasmid_dir()` is searched second — it extends the library rather than overriding it. This means pYTK part names are stable regardless of what is in your personal collection, and you can point `set_plasmid_dir()` at your assembly output folder to make newly-built plasmids immediately available for the next round of assembly without copying anything.

To add a plasmid permanently to the bundled library, drop its `.gb` or `.fa` file into `inst/extdata/plasmids/` and reinstall the package.

## Batch script: `build_entry_plasmids.R`

A driver script for the common case of building N Type 3 entry plasmids from N PCR products, all going into the same destination vector. Lives at the package root. Depends on the installed package — no path configuration needed.

### Input format

A tab-separated file (default name `gg_parts_summary.tsv`) with these required columns:

| Column | Type | Notes |
|---|---|---|
| `gene_name` | string | Used to name the output file (`pYTK_T3_<gene_name>.gb`); non-filesystem-safe characters are replaced with `_` |
| `species` | string | Carried into the summary TSV; not otherwise used |
| `accession` | string | Carried into the summary TSV; not otherwise used |
| `fwd_primer` | string | Carried; not otherwise used |
| `rev_primer` | string | Carried; not otherwise used |
| `codon_optimized_sequence` | string | Carried; not otherwise used |
| `pcr_product` | string | The full top-strand sequence of the linear PCR product, including primer tails. Must contain BsaI recognition sites positioned so that digestion releases a fragment with overhangs matching the destination backbone. |

The script reads `pcr_product` as the linear amplicon, digests it with BsaI, picks the largest released fragment as the insert, and ligates it into a freshly-digested backbone fragment. The other six columns are passed through untouched into the summary output.

### Configuration

Two variables at the top of the script:

```r
DESTINATION_PLASMID <- "pYTK_T3_dest_demo"   # plasmid name, resolved via load_plasmid()
USER_PLASMID_DIR    <- NULL                   # optional supplemental plasmid directory
ENZYME              <- "BsaI"
```

`DESTINATION_PLASMID` is a name (not a path); it's resolved through the library search path. Set `USER_PLASMID_DIR` to point at a folder of custom plasmids — the bundled pYTK library is always searched first.

### Run

```bash
# Default: reads gg_parts_summary.tsv, writes to entry_plasmids/
Rscript build_entry_plasmids.R

# Explicit paths:
Rscript build_entry_plasmids.R my_parts.tsv my_outputs/
```

### Output

The script creates `<output_dir>/` (default `entry_plasmids/`) containing:

- One GenBank file per successful row, named `pYTK_T3_<gene_name>.gb`
- `ligation_summary.tsv` with the original input columns plus four new ones:

| Column | Values |
|---|---|
| `status` | `OK` (clean ligation), `OK_WITH_WARNING` (succeeded but the package noted something — typically internal recognition sites in the product, which is normal for Type 3 entry plasmids), or `FAIL` (assembly impossible) |
| `product_bp` | Length of the assembled circular product, or empty for `FAIL` rows |
| `output_file` | Path to the GenBank file, or empty for `FAIL` rows |
| `message` | One-line description: overhang summary for OK, warning text for OK_WITH_WARNING, error message for FAIL |

The script exits with status 1 if any row failed, so it's safe to chain in a Makefile or shell pipeline. Per-row failures don't abort the batch — every row is attempted.

### Design assumption

`select = "largest"` is hardcoded for both the destination digestion and the PCR-product digestion. This is correct when:

- The destination's largest fragment is the CamR + ori backbone (true for any sensible YTK Type 3 destination — the GFP dropout is much smaller).
- The PCR-product's largest fragment is the released CDS (true when the primer tails follow the standard `CGAGCG-GGTCTC-N-[overhang]-[anneal]` design — the primer-tail edge fragments are about 11–15 bp each, much smaller than any real CDS).

If your primer tails are unusually long or your gene is under ~50 bp, edit the script to use `select = c(min_bp, max_bp)` instead.

## Batch assembly: `assemble_interactive` and `assemble_from_table`

Two convenience functions for running multi-part cassette assemblies without writing boilerplate `load_plasmid` + `golden_gate_assemble` calls.

### Interactive

Prompts for a plasmid name and a space-separated part list, then assembles and writes the GenBank file:

```r
assemble_interactive(enzyme = "BsaI", output_dir = "./plasmids")
# Output plasmid name: MY_CONSTRUCT
# Parts (space-separated): pYTK090 pYTK008 pYTK047 pYTK073 pYTK074 pYTK086 pYTK092
```

### Table-driven

Reads a TSV or CSV manifest — one row per assembly, one column per part slot — and writes a GenBank file for each row. An example manifest (`assemblies_example.tsv`) is included at the package root.

Required column: `name` (output plasmid name).

Optional part columns (empty cells are skipped):

`vec` | `1` | `2` | `3` | `3a` | `3b` | `234` | `4` | `4a` | `4b` | `5` | `6` | `7` | `8` | `8a` | `8b` | `678`

```r
results <- assemble_from_table("assemblies_example.tsv", output_dir = "./plasmids")
```

Returns a data frame with `name`, `status` (`OK` / `FAILED` / `SKIPPED`), `output`, and `error` columns. Failures are reported per-row and do not abort the batch.

To use newly-assembled plasmids as inputs for a subsequent round of assembly, point `set_plasmid_dir()` at the same output folder:

```r
set_plasmid_dir("./plasmids")
results <- assemble_from_table("level2_assemblies.tsv", output_dir = "./plasmids")
```

## Testing

The package ships with five test files at the package root:

```bash
Rscript test_run.R              # 9 unit tests on synthetic plasmids
Rscript test_features.R         # End-to-end: GenBank read -> assemble -> GenBank write
Rscript test_ligate.R           # 9 tests on the gel-purification + ligation workflow
Rscript test_origin_wrap.R      # Regression test for origin-spanning cuts
Rscript test_plasmid_library.R  # 7 tests on load_plasmid / list_plasmids / set_plasmid_dir
```

## Limitations

- Type IIS only — no support for blunt or sticky Type II enzymes (EcoRI, HindIII, etc.).
- Returns the unique circular product or errors. Does not enumerate linear or partial products.
- Features with `join(...)` or `order(...)` locations in input GenBank files are dropped.
- The GenBank parser is hand-rolled; it handles SnapGene/Benchling/NCBI files in testing but exotic edge cases may need work.

## Comparison with `moclo-ytk`

The Python `moclo-ytk` library does considerably more than this package: it validates that each part conforms to the YTK Type 1–8 schema (correct overhangs, expected internal features), maintains a registry of curated Addgene parts, and supports the full Level 0 → Level 1 → Level 2 tier hierarchy. `goldengateR` is deliberately narrower: a digestion-and-ligation simulator with a small plasmid library and a batch driver. If you need YTK-specific schema validation, use `moclo-ytk` for the design-time check; use `goldengateR` for the in-pipeline simulation step where staying in R matters (e.g., generating reference sequences for downstream Nanopore read alignment).
