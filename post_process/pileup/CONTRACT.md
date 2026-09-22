# Pileup provider contract, version 1

A **provider** turns (alignments, sites) into per-site read counts. An **annotator** puts those counts
into a VCF. This file is the only thing they share.

The split exists so the pileup engine can be replaced — an accelerated C++ implementation is expected
— without touching the annotator, and so a replacement can be proven equivalent before it is trusted
(see `tests/pileup/`).

## Scope

A provider is given one sample's alignments and a list of sites, and reports what the reads say at
each site. It does **not** decide anything: no filtering, no genotyping, no thresholds. Interpretation
belongs to the annotator and to the filters above it.

## Invocation

```
run_pileup.sh --provider NAME --alleles FILE --regions FILE --ref FASTA --bam FILE \
              --role ROLE [--depth N] [--min-mq N] [--min-bq N] [--baq on|off] > canonical.tsv
```

One invocation handles **one sample**. This is deliberate and is not an efficiency oversight: the
bias statistics `bcftools` computes (`BQBZ`, `MQBZ`, `MQSBZ`, `RPBZ`, `SCBZ`, `MQ0F`, `SGB`, `VDB`,
`FS`) are INFO-level and pool every sample in the call, so a multi-sample pileup makes subject-
specific bias unrecoverable. Per-sample invocation also removes any dependence on sample naming or
column order — the role comes from `--role`, never from a parsed `@RG SM` field.

## Output format

A tab-separated table on stdout, sorted by `CHROM` then `POS`, with a preamble of `#`-prefixed lines:

```
#contract_version=1
#provider=<name>:<version>
#params=depth=10000,minMQ=0,minBQ=13,baq=off
#role<TAB>file<TAB>kind<TAB>vcf_sample_name
```

The `#role` line records what this invocation covered, so the annotator never has to infer the
mapping from sample names. `kind` is `bam`, `cram` or `altonly` (the `Mutations.fastq.bam` fallback,
which contains no reference reads).

Then one header line naming the columns, then one row per site:

| # | column | meaning |
|---|---|---|
| 1 | `CHROM` | |
| 2 | `POS` | VCF 1-based, on the **normalized** coordinates the caller supplied |
| 3 | `REF` | exactly as supplied |
| 4 | `ALT` | exactly as supplied |
| 5 | `ROLE` | `SUBJECT`, `CONTROL_1`, … |
| 6 | `DP` | reads considered at the site |
| 7-8 | `AD_REF` `AD_ALT` | reads supporting each allele |
| 9-12 | `ADF_REF` `ADF_ALT` `ADR_REF` `ADR_ALT` | the strand 2x2 |
| 13-16 | `F1R2_REF` `F1R2_ALT` `F2R1_REF` `F2R1_ALT` | read orientation; `.` in v1 |
| 17 | `MQ0F` | fraction of reads with mapping quality 0 |
| 18-23 | `RPBZ` `BQBZ` `MQBZ` `MQSBZ` `SCBZ` `SGB` | bias statistics |
| 24 | `SP` | phred-scaled strand-bias p-value |
| 25 | `SCR` | soft-clipped reads at the site |
| 26 | `NMBZ` | mismatch-count bias |

## Rules a provider must satisfy

These are what `tests/pileup/` checks, and what makes a replacement safe.

1. **One row per supplied site, always.** A site with no read support is reported with zeros, not
   omitted. A missing row and a zero row mean different things — "the provider did not look" versus
   "the reads say nothing" — and only the second is evidence. With `bcftools` this requires
   `bcftools call -i`; without it, unsupported sites vanish silently.

2. **`REF` and `ALT` come back exactly as supplied.** Counts must be attributable to the caller's
   allele, not to whatever the engine decided the site contains. `bcftools mpileup` rewrites indels
   into full repeat context — `A>AT` becomes `ATTTTTT>ATTTTTTT` — which broke 18 of 49 indels in
   testing before a trailing `bcftools norm` restored the input representation. A provider owns this
   normalization; the annotator will not re-derive it.

3. **`AD_REF`/`AD_ALT` refer to the supplied REF and ALT, in that order.** Not to an engine-chosen
   allele ordering.

4. **Unsupported fields are `.`, never a guessed value, and never a fatal error.** A provider that
   cannot compute soft-clip bias emits `.` for it. The annotator degrades; it does not fail. The v1
   column set deliberately includes `F1R2`/`F2R1`, which the bcftools provider cannot produce and
   emits as `.` — so this path is exercised by the only provider that exists, rather than being
   untested until the day it matters.

5. **No interpretation.** No filtering on depth, no dropping of low-quality sites, no genotype calls.

6. **Deterministic.** The same inputs produce the same bytes, modulo the `#provider` and `#params`
   preamble lines.

## `*` alleles

Sites whose ALT is the spanning-deletion placeholder `*` are not piled up: `*` denotes absence of
sequence, so there is nothing to count. The caller omits them from `--alleles`, and the annotator
emits their row with a reason code rather than silently dropping them — the VCF should say that no
evidence was sought, not imply that none was found.

## Known limits of the bcftools provider

A provider reports what its engine can score, and `bcftools mpileup` cannot score every allele RUFUS
can call. Two classes come back with `AD_ALT=0` **even though the reads plainly carry the variant**:

| class | example from `tests/pileup/fixture` | what the reads show |
|---|---|---|
| composite / MNV | `chr20:150000 ACA>CAC` | 13 of 32 reads carry all three substitutions; scored 0 |
| large deletion | `chr20:170000`, 1000 bp | 127 reads span the deleted span; scored 0 |

mpileup has no MNP model, so a composite allele cannot be counted as one unit; and a deletion far
longer than a read is not represented in a per-position pileup at all.

**This is the single most important thing to know before filtering on these numbers.** For these
classes `AD_ALT=0` means "this engine cannot answer", NOT "there is no support" — the two are
indistinguishable in the output, and treating the first as the second would reject exactly the
variants an assembly-based caller exists to find. Consumers must gate on variant class, never on
`AD_ALT` alone. SNVs and short indels are scored correctly and carry no such caveat.

The annotator marks the detectable half. A **composite allele is structural** — REF and ALT the same
length and longer than one base — so `INFO/NO_PILEUP_MODEL=1` is set on those with no threshold and no
guessing, and `FORMAT/AD_OTHER` reports reads supporting neither listed allele (14 of 33 at the
fixture's MNV: the variant-carrying reads themselves).

**Large deletions are NOT flagged**, and this is a known gap rather than an oversight. Detecting them
needs a deletion-length-versus-read-length threshold, which is a judgement call rather than a
structural fact — and unlike the MNV they give no secondary signal either: at the fixture's 1000bp
deletion `AD_OTHER` is 0, because the spanning reads align cleanly as reference. So a long deletion
still reads as `AF=0, NO_PILEUP_MODEL=0`, indistinguishable from a variant that is genuinely absent.
Until that threshold is chosen, gate on variant class.

The conformance fixture deliberately contains one of each, asserted as blind. If a future provider
scores them, the test fails — which is the point: an improvement should be noticed and the contract
updated, not absorbed silently.

## Versioning

`#contract_version` changes when a column's meaning changes or a column is removed. Adding a column
at the end, or filling in one that was previously always `.`, does not require a bump, because rule 4
means consumers already tolerate `.`.

## Conformance

`tests/pileup/` holds a small fixture and a golden table. A provider is conformant when it reproduces
the golden and satisfies the rules above. `tests/pileup/compare_providers.sh` diffs two providers over
the same input with a numeric tolerance, which is how a replacement is qualified before it is trusted
with real data.
