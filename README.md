# iedb-epitope-mapping-nextflow

Nextflow pipeline that maps immune epitopes from the [Immune Epitope Database (IEDB)](https://www.iedb.org/) onto a reference proteome.

## Overview

IEDB catalogs experimentally characterized immune epitopes (short peptide sequences recognized by antibodies or T cells) and their source proteins, but a given epitope's exact position on a particular organism's reference proteome isn't provided directly. This pipeline takes an IEDB epitope table and a reference proteome for a target organism/taxon, determines which epitopes originate from proteins in that taxon (using NCBI Taxonomy and GenBank), matches epitope peptide sequences against the reference proteome using [PEPMatch](https://github.com/IEDB/PEPMatch), and produces an indexed GFF3 of epitope-to-protein mappings for display in the VEuPathDB genome browser.

Peptides are split into three categories, each searched with parameters suited to it:
- very short peptides (matched with a smaller k-mer size, since short peptides need finer-grained indexing),
- peptides whose source protein belongs to the reference taxon (matched allowing one mismatch, since these are expected to align closely but not necessarily perfectly), and
- all other peptides (matched exactly).

Peptides whose source taxon matches the reference are additionally checked for an exact full-length match against their original GenBank source protein.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- A container engine — Docker or Singularity. `nextflow.config` includes `conf/docker.config` by default (`docker.enabled = true`); `conf/singularity.config` and `conf/lsf.config` (Singularity plus LSF cluster execution) are available as alternative profiles.

## Usage

```
nextflow run VEuPathDB/iedb-epitope-mapping-nextflow -r main \
  --refFasta /path/to/reference_proteome.fsa \
  --peptidesTab /path/to/iedb_epitopes.tab \
  --taxon 5833 \
  --results /path/to/output \
  -resume -C my.config
```

The pipeline has a single, unnamed entry point (`workflow { ... }` in `main.nf`, which calls the `epitopeMapping` workflow), so no `-entry` flag is needed.

Steps performed:
1. `fetchTaxon` retrieves every descendant taxon ID under `params.taxon` from NCBI Taxonomy (via EDirect).
2. `peptideProteinAccessionsFilteredByTaxa` filters the epitope table against that taxon list, splitting peptides into taxon-matching and non-matching sets and separating out very short peptides, while also collecting the accessions of their source proteins.
3. `fetchProtein` downloads the full source protein sequences for those accessions from GenBank (via EDirect).
4. `iedbExactMatches` checks taxon-matching peptides for an exact full-length match against their source protein and the reference proteome.
5. Three parallel invocations of the `pepMatch` subworkflow (`smallExactPepMatch`, `exactPepMatch`, `inexactForTaxaPeptidesPepMatch`) each preprocess the reference proteome with PEPMatch (`pepmatch-preprocess`) and match their respective peptide subset against it (`pepmatch-match`), filtering the results down to the relevant columns.
6. `mergeResultsFiles` combines the PEPMatch results and the exact-match results into a single GFF, using `params.nonTaxaShortPeptideCutoff` to decide how short non-taxon-matching peptide matches are treated.
7. `indexResults` sorts the merged GFF, compresses it with `bgzip`, and indexes it with `tabix`, publishing the result to `params.results`.

## Key Parameters

| Parameter | Description | Default |
|---|---|---|
| `params.refFasta` | Path to the reference proteome FASTA for the target organism | `data/AnnotatedProteins.fsa` |
| `params.peptidesTab` | Path to the IEDB epitope tab file | `data/iedb.tab` |
| `params.taxon` | NCBI taxonomy ID of the reference organism | `5833` |
| `params.peptideMatchResults` | Filename for the merged GFF prior to indexing | `iedb.gff` |
| `params.results` | Directory the final indexed GFF is published to | `output` (relative to launch directory) |
| `params.nonTaxaShortPeptideCutoff` | Minimum length for a non-taxon-matching short peptide to be kept in the merged results | `5` |

Per-subworkflow PEPMatch parameters (k-mer size, index format, allowed mismatches) are set via `process.withName` selectors in `nextflow.config` for each of the three peptide categories and are not typically changed per run.

## Output

A sorted, `bgzip`-compressed and `tabix`-indexed GFF3 file (named per `params.peptideMatchResults`, e.g. `iedb.gff.gz` plus its `.tbi` index), published to `params.results`. Each feature represents an IEDB epitope mapped onto the reference proteome, annotated with whether it was an exact or approximate match and whether its source protein belongs to the reference taxon.
