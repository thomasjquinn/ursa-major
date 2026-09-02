# Commit 13: Bundled Example Data

## Commit Abstract: Carry the fifteen example data files bundled with upstream baerhunter into this fork unchanged, so that the package's man pages and vignette have something to build against.

### Commit Number:

`4e19f08307f1055c255ba56040447c871208b5e4`, short form `4e19f08`. Parent `78c263e`. Landed 31 August 2026.

### Commit Summary:
Add upstream example data for package examples

### Commit Description:

```
data(inst): carry upstream inst/extdata unchanged

Add the fifteen files bundled with upstream baerhunter: six
BAM subsets covering positions 1 to 10,000, their indexes, the
EnsemblBacteria release 40 H37Rv annotation, and two small
text files.

Carried whole and unmodified from irilenia/baerhunter at
1e64171. They supply the example data the package's man pages
and vignette build against, and were absent from this fork
while it was maintained as a set of source files rather than
an installable package.

Source hashes recorded before the copy and verified against a
fresh clone afterwards; samtools quickcheck passes on all six
alignments.
```

### Commit Details:

Fifteen files added under `inst/extdata/`, 20,610 insertions and no deletions. Nothing already in the repository is changed.

**No source file is touched and no code is written.** This commit adds data only, carried whole from `irilenia/baerhunter` at `1e64171`. It is recorded as its own commit because the files are large, binary, and reviewed by provenance rather than by reading.

---

## What was added

**Issue Summary:** This fork was maintained as a contribution overlay holding only the source files under development, so upstream's `inst/extdata` was never carried across. A package whose man pages and vignette build against example data cannot install and check without it.

**Solution Summary:** Copy the fifteen files whole from upstream and add them under `inst/extdata/`, with no renaming, no subsetting and no re-indexing.

**Note:** This is a completeness fix rather than a change of any kind. Every file is byte-identical to its upstream original, and nothing in the package reads them until the man pages and the vignette exist.

### The file tree

```
inst/extdata/
├── ERR262980_1_10000.bam                          163,609 bytes
├── ERR262980_1_10000.bam.bai                           96 bytes
├── ERR262982_1_10000.bam                          170,387 bytes
├── ERR262982_1_10000.bam.bai                           96 bytes
├── ERR262983_1_10000.bam                          180,561 bytes
├── ERR262983_1_10000.bam.bai                           96 bytes
├── ERR262984_1_10000.bam                           63,764 bytes
├── ERR262984_1_10000.bam.bai                           96 bytes
├── ERR262987_1_10000.bam                           79,497 bytes
├── ERR262987_1_10000.bam.bai                           96 bytes
├── ERR262988_1_10000.bam                           76,172 bytes
├── ERR262988_1_10000.bam.bai                           96 bytes
├── Mycobacterium_tuberculosis_h37rv.ASM19595v2.40.chromosome.Chromosome.gff3
│                                                2,395,583 bytes
├── chromosome.txt                                      22 bytes
└── conditions.txt                                     143 bytes
```

**Six BAM subsets**, each covering positions 1 to 10,000 of the H37Rv chromosome, from the Cortes et al. RNA-seq accessions `ERR262980`, `ERR262982`, `ERR262983`, `ERR262984`, `ERR262987` and `ERR262988`. Three are exponential-phase and three stationary-phase, as `conditions.txt` records.

**Six BAI indexes**, one per alignment, uniformly 96 bytes. The uniform size is expected rather than suspicious: each index covers a 10,000-base subset of a single chromosome, so the linear index is minimal and the same shape in every file. This was confirmed rather than assumed, below.

**One GFF3 annotation**, the EnsemblBacteria release 40 H37Rv chromosome annotation, at 2.4 MB the largest file in the commit and 96% of its total size.

**Two small text files.** `chromosome.txt` is a one-line alias mapping the annotation's sequence name to the accession, `Chromosome,AL123456.3`. `conditions.txt` is a seven-line sample sheet giving each BAM its condition label.

---

## Verification

**Issue Summary:** Binary files cannot be reviewed by reading, and the preceding commit's `.gitattributes` rules had never been exercised, so a line-ending-normalised BAM would have been invisible in the diff and would have produced different coverage without erroring.

**Solution Summary:** Record the source hashes before the copy, confirm the binary rules took effect at staging, run `samtools quickcheck` on all six alignments, and verify the round trip by re-cloning and checking every file against its recorded hash.

**Note:** The round trip is the condition that matters. It proves that the files reached GitHub and came back byte-identical, which is the only evidence that the previous commit's ordering did what it was written to do.

### Conditions and results

| condition | result |
|---|---|
| fifteen files staged, not fourteen | 15 |
| nothing excluded by `.gitignore` | confirmed; no exception lines were needed |
| BAMs treated as binary | GitHub Desktop reported *"This binary file has changed"* rather than rendering them as text |
| `samtools quickcheck` on all six alignments | passed, samtools 1.23.1 |
| round trip against the recorded hashes | `sha256sum -c` returned `OK` fifteen times from a fresh clone |

**The binary-rule confirmation is the observable one.** Had `text=auto` won, the alignment files would have been rendered as text at staging and normalised on checkout. They were not, which is the evidence that the rules committed in the previous commit were in force by the time they were needed.

**`samtools quickcheck` also settled a standing question** about the six uniform 96-byte indexes: all six passed, so the uniformity is a property of the data rather than a sign of truncated or corrupt index files.

### Recorded hashes

sha256, first sixteen characters, in tree order.

```
d5852cdcdf89c3f5  ERR262980_1_10000.bam
1c7ce26b3343d9ef  ERR262980_1_10000.bam.bai
51725c613966ab6d  ERR262982_1_10000.bam
9e3397ac6abfc560  ERR262982_1_10000.bam.bai
8060d087c4339a27  ERR262983_1_10000.bam
742848abb621abe6  ERR262983_1_10000.bam.bai
f8abec861811509b  ERR262984_1_10000.bam
7888be47667036a6  ERR262984_1_10000.bam.bai
2a85c8a6e6acbd5d  ERR262987_1_10000.bam
17357a1b123ecf53  ERR262987_1_10000.bam.bai
935a5fb82c91ce6b  ERR262988_1_10000.bam
ea42de7573d9d6f8  ERR262988_1_10000.bam.bai
9a313f2940d12a60  Mycobacterium_tuberculosis_h37rv.ASM19595v2.40.chromosome.Chromosome.gff3
954adc2534c508ec  chromosome.txt
0512db19fcbf6a2b  conditions.txt
```
