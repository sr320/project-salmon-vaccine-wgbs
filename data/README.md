# data

Raw whole-genome bisulfite sequencing (WGBS) reads from Atlantic salmon (*Salmo salar*) parr sampled before and after vaccination, plus a sample metadata table.

> **Draft.** Items marked **TODO** could not be determined from the files themselves and need to be filled in by someone who knows the experiment.

## Contents

| File | Description |
|---|---|
| `Metadata.tsv` | Sample sheet: library, read, condition, timepoint, replicate (one row per FASTQ file) |
| `SsParr*_1.fastq.gz` | Read 1 of each paired-end library |
| `SsParr*_2.fastq.gz` | Read 2 of each paired-end library |

Total size is about 31 GB. The FASTQ files are too large for GitHub and should not be committed (they are excluded by the repository `.gitignore`).

## Experimental design

Four libraries: two conditions with two biological replicates each.

| Library | Condition | Timepoint | Replicate |
|---|---|---|---|
| `SsParrT0_1` | Before vaccination | 0 ATUs | 1 |
| `SsParrT0_2` | Before vaccination | 0 ATUs | 2 |
| `SsParrVac_1` | After vaccination | 600 ATUs | 1 |
| `SsParrVac_2` | After vaccination | 600 ATUs | 2 |

ATU = accumulated thermal units (degree-days) since vaccination.

**TODO:** tissue sampled, vaccine used, whether each library is one fish or a pool, whether the same fish were sampled at both timepoints, fish source and rearing conditions, sampling dates.

## File naming

```
SsParr{T0|Vac}_{replicate}_{read}.fastq.gz
```

- `Ss` = *Salmo salar*, `Parr` = life stage
- `T0` = before vaccination, `Vac` = after vaccination (600 ATUs)
- first number = biological replicate (1 or 2)
- second number = read of the pair (1 = R1, 2 = R2)

So `SsParrVac_2_1.fastq.gz` is R1 of post-vaccination replicate 2.

## Sequencing files

Paired-end, 2 × 151 bp, Illumina (instrument `LH00714`; quality scores are binned). Read counts are read pairs; R1 and R2 counts match for every library.

| Library | Read pairs | Flowcell | Lane | Index (i7+i5) |
|---|---:|---|---:|---|
| `SsParrT0_1` | 77,603,127 | 22MGM5LT4 (run 97) | 1 | TGTCGGAT+TGGTGAGC |
| `SsParrT0_2` | 68,776,693 | 22MG7GLT4 (run 94) | 6 | CGCTGGTC+ATACGATA |
| `SsParrVac_1` | 69,143,094 | 22MG7GLT4 (run 94) | 6 | ATTGGATA+AGGTAGCC |
| `SsParrVac_2` | 70,172,391 | 22MG7GLT4 (run 94) | 6 | GATGGTAC+CTTGAACT |

Flowcell, lane, and index are the values seen in the headers of the first 100,000 reads of each R1 file (one flowcell and lane per file; the index shown is the dominant one, at 98% or more of reads).

Each library is roughly 21 to 23 Gb of sequence, about 7× to 8× raw coverage of the ~3 Gb Atlantic salmon genome before trimming, deduplication, and mapping losses.

**TODO:** library prep kit and conversion method, sequencing facility, where the raw data are archived (e.g. SRA/ENA accession or lab storage URL), any spike-in controls (lambda, pUC19) for estimating conversion efficiency.

### Checksums (MD5)

```
2ccbc4fdf95ad060948cbf72d856a5f9  SsParrT0_1_1.fastq.gz
09bb5ebf63d6735d51e18d5806f2c135  SsParrT0_1_2.fastq.gz
636630d2825135867e61f7669716773a  SsParrT0_2_1.fastq.gz
6075c2f59477ba7ab88e096c2060fad3  SsParrT0_2_2.fastq.gz
4b76167847b3814c7354208f555a2891  SsParrVac_1_1.fastq.gz
7875b9e02461f016458e68d80752236e  SsParrVac_1_2.fastq.gz
57da1e5e2e4a55899319d1053ffca602  SsParrVac_2_1.fastq.gz
88ff9496b3b39b3af584ed7e661e61d9  SsParrVac_2_2.fastq.gz
```

These were computed on the local copies on 2026-11-25. They have not been compared against checksums from the sequencing provider.

## Metadata.tsv

Tab-delimited, 8 data rows, one per FASTQ file. Group by `Library` to get the 4 samples.

| Column | Description |
|---|---|
| `SampleID` | FASTQ file name without `.fastq.gz` (includes the read number) |
| `Library` | Library the file belongs to; R1 and R2 of a pair share this value |
| `Read` | `R1` or `R2` |
| `Condition` | `Before vaccination` or `After vaccination` |
| `Timepoint` | `0 ATUs` or `600 ATUs` |
| `Replicate` | Biological replicate, 1 or 2 |

The file was cleaned on 2026-11-25: Windows line endings and a leading space in the `After vaccination` values were removed, a final newline was added, `Library` and `Read` columns were added, and `replicate` was renamed `Replicate`. The original had one row per file with nothing linking the two reads of a pair.

## Notes

From a review of the files on 2026-11-25:

1. **Batch is partly confounded with condition.** `SsParrT0_1` was sequenced on a different flowcell and run than the other three libraries. With two replicates per group, a run effect on that library cannot be separated from biological variation.
2. **Reads are directional and converted as expected.** In the first 100,000 reads of each file, cytosine makes up about 2% of R1 bases and about 21% of R2 bases, consistent with a directional bisulfite-type library.
3. **R2 has a 5′ artifact.** R2 reads frequently begin with a `GGGAG`/`GGGA`-like motif, and at least some carry adapter sequence shortly after. Check the per-base content plots in FastQC and expect to hard-trim the 5′ end of R2 (and likely the 3′ end of R1) in addition to adapter trimming. The right number of bases depends on the library kit (TODO above).
4. **Small design.** n = 2 per group limits differential methylation analysis to large, consistent effects.
5. **FASTQ files are git-ignored.** The repository `.gitignore` excludes `data/*.fastq.gz`; only this README and `Metadata.tsv` are meant to be tracked.
