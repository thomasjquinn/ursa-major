# Parameter Scout: Instructions

For Parameter Scout **Version 1.0**

---

## **What is Parameter Scout?**

Parameter Scout is a simple R utility that helps you choose coverage cutoff parameters. This is an optional step you may run before the main baerhunter pipeline.

Picking the correct cutoff parameters often seems more like an art than a science. Parameter Scout measures the intergenic coverage in the BAM files and suggests two pairs of cutoff parameters based on the data itself.

baerhunter (Ozuna et al., 2020) requires both `low_coverage_cutoff` and `high_coverage_cutoff` to determine which parts of your data become predicted features and which are ignored. These are raw read depths, unscaled to your library or genome, so a cutoff that works well with one dataset can be a poor choice for another. For example, a low coverage cutoff of 5 and high coverage cutoff of 10 might be sensible for a shallow library, but far too low for a deep one.

Parameter Scout is only an advisory tool. It does not write an annotation file. It does not automatically set any parameter within baerhunter. Instead, it provides a picture of the data and possible coverage cutoff parameters. It does not pass anything to the rest of the pipeline. You still need to choose the cutoff parameters and edit `feature_file_editor()` yourself.

If you already know what `low_coverage_cutoff` and `high_coverage_cutoff` your data needs, you do not need to run Parameter Scout at all.  

### There are two versions, and you need the right one

Parameter Scout ships as two files:



| file | use it when | the core function |
| --- | --- | --- |
| `parameter_scout_paired_end.R` | your reads are paired-end | `parameter_scout_paired_end()` |
| `parameter_scout_single_end.R` | your reads are single-end | `parameter_scout_single_end()` |



The two files are **almost identical**. They differ only in fifteen lines: the header, the name of the function you call, and one default setting. Everything else is the same code.

---

## **How Parameter Scout Works**

### What it measures

Parameter Scout reads each of your BAM files once. For each one, it works out the coverage at every base, then throws away everything except the **intergenic regions**, the parts of the genome between the annotated features.

The intergenic regions are worked out from your GFF file, using the same rules the main pipeline uses. In Parameter Scout, tRNA and rRNA genes are treated as features, so they are **excluded** from the background. This matters, because rRNA coverage can be enormous and would drag every number upwards.

### What the numbers describe

Parameter Scout measures coverage only in the intergenic regions, the parts of the genome between the annotated features. Within those regions it looks only at the positions where at least one read landed. Positions with no reads are left out.

This matters when you read the median. Most of the intergenic space in a bacterial genome carries no reads at all. If those empty positions were counted, the median would sit at or near zero for almost any dataset, and it would tell you nothing useful about where noise ends. Leaving them out means the median answers a more useful question: among the intergenic positions that show any signal, what is a typical amount? So a median of 15 means that half of the covered intergenic positions have 15 reads or more. It does not mean that half of the intergenic genome is covered at that level.

### From percentiles to a suggestion

For each BAM file, Parameter Scout reports six percentiles: P25, the median, P75, P80, P85 and P90. You get one row per file.

To turn many rows into one number, Parameter Scout takes the **median across files** of each column, then rounds up to the next whole number. Rounding up rather than to nearest is deliberate.

From that aggregate row, two pairs are offered:



| pair | low cutoff | high cutoff |
| --- | --- | --- |
| **Relaxed** | Median | P75 |
| **Stringent** | P75 | P80 |



The two sets of parameters (Median & P75; P75 & P80) follow the method outlined in Sivasankaran (2024).  
Notice that **P75 appears twice**: it is the high cutoff of the Relaxed pair and the low cutoff of the Stringent pair. One number is doing two jobs. Do not use it for both: a low and a high cutoff set to the same value is not one of the pairs offered here, and not a valid choice.

### The two cutoffs do not behave the same way

This is the single most useful thing to understand before you choose your parameter pair.

**The low cutoff decides where features start and stop.** The coverage vector is sliced wherever it crosses this level, and each slice becomes a candidate feature. Raise the low cutoff and the slice boundaries move: features get shorter, and one long feature can break into several shorter ones, each of which may then pass. So raising the low cutoff **can give you more features, not fewer.**

**The high cutoff decides how much evidence a feature needs.** It does not move any boundary. It only asks whether enough of the slice sits above the level. Raise it and you get a strict subset of what you had before (the same features, minus some).

Only the second behaves the way you would expect. **A reader who assumes "Stringent" means fewer predictions will sometimes be wrong**, because the Stringent pair raises the low cutoff as well, and that part is not predictable in the same way.

Parameter Scout is a starting point for exploring your data. Other cutoffs may suit your data better than what Parameter Scout suggests.

---

## **How to Use Parameter Scout**

### What you need installed first

Parameter Scout **attaches** these packages for you when you load it. It cannot **install**
them. If any one is missing, loading the module stops with a message naming it, so it is
worth checking before you start.



| package | where it comes from |
| --- | --- |
| `S4Vectors` | Bioconductor |
| `IRanges` | Bioconductor |
| `GenomicRanges` | Bioconductor |
| `GenomicAlignments` | Bioconductor |
| `Rsamtools` | Bioconductor |
| `Seqinfo` | Bioconductor. On older Bioconductor releases these functions live in `GenomeInfoDb` instead, and Parameter Scout attaches whichever you have |
| `ggplot2` | CRAN |
| `grid` | already installed: it comes with R |



To check what you already have, run this in the RStudio Console:



```
need <- c("S4Vectors", "IRanges", "GenomicRanges", "GenomicAlignments",
          "Rsamtools", "ggplot2")
need <- c(need,
          if (requireNamespace("Seqinfo", quietly = TRUE)) "Seqinfo"
          else "GenomeInfoDb")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
cat("R", format(getRversion()), "\n")
if (length(missing)) cat("MISSING:", paste(missing, collapse = ", "), "\n") else
  cat("All packages found.\n")
```



To install the ones you are missing:



```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("S4Vectors", "IRanges", "GenomicRanges",
                       "GenomicAlignments", "Rsamtools", "Seqinfo"))

install.packages("ggplot2")
```



If `BiocManager` cannot find `Seqinfo`, your Bioconductor is an older one and the same
functions are in `GenomeInfoDb`. Install that instead. Parameter Scout checks which of the
two you have and uses it, so you do not need to tell it which.

Parameter Scout was built and tested on R 4.6.0 with Bioconductor 3.23, which is
the setup to use if you have the choice. Its test suite also passes on R 4.3.3, so an
older R may well work; no lower bound has been established.

### Step 1\. Put everything in one folder, then point R at that folder

Use one folder for one dataset. It should hold:

* the BAM files for that dataset, and no others
* the GFF annotation file
* `count_features.R`, `feature_file_editor.R`, and the Parameter Scout file (parameter\_scout\_paired\_end.R or parameter\_scout\_single\_end.R)

Then make that folder your working directory:



```
setwd("C:/Users/you/Documents/mydata")
getwd()    # check it worked
```



In RStudio you can do the same from the menu: **Session > Set Working Directory > Choose Directory**.

**That one line does three jobs, and skipping it causes three different problems.**

1. **It lets step 2 find the R files.** `source("count_features.R")` looks for the file in the working directory. If R is looking somewhere else, the line fails and nothing can run.
2. **It decides where Parameter Scout writes.** By default the results go into a folder called `scout` inside your working directory. Parameter Scout creates that folder for you, so there is nothing to make first.
3. **It decides where baerhunter writes.** baerhunter never infers its output location from where your BAM files are. Every output filename you give it is resolved against the working directory. Setting it once, here, means `feature_file_editor()` and `count_features()` write beside your data instead of somewhere you did not expect. baerhunter also does not create folders, so any folder you name in an output path has to exist already.

**Parameter Scout reads every `.bam` file in the folder, whatever the case of the extension.** It does not read a list you have in your head. If a seventh BAM appears in a six-BAM folder, your percentiles change and nothing warns you.

The same is true of `feature_file_editor()`, which finds its BAM files the same way. So the folder _is_ a parameter, and keeping it tidy is part of getting the right answer. If you want to scan a subset, use the `bam_txt_list` argument to name the files instead of moving them about.

The run log records exactly which files were read, so you can check afterwards.

**How many at once?** The tool is designed for **1 to 9 BAM files**. It stops with a message above nine.

* If you want to **read the shape** of the distributions, scan **no more than 6 BAM files**. The figure is a fixed height whatever the number of files, so each distribution is drawn shorter as the count goes up. At nine they are too squashed to read, and reading the shape is what the figure asks you to do.  
* If you only want the **numbers**, nine BAM files is fine. The table and the suggested pairs are text, so nothing is lost.

### Step 2\. Load the three files

Three lines. There are no `library()` calls here because Parameter Scout attaches the
packages listed above itself, as soon as you load it.



```
source("count_features.R")           # supplies .resolve_gff_cache()
source("feature_file_editor.R")      # supplies major_features()
source("parameter_scout_paired_end.R")
```



All three files must be in your working directory, which is what step 1 set. Use
`parameter_scout_single_end.R` on the third line if your reads are single-end.



### Step 3\. Run the scan

Only two arguments have to be supplied. Everything else has a default.

**Call the function that matches the file you loaded in step 2**:
`parameter_scout_paired_end()` for paired-end reads,
`parameter_scout_single_end()` for single-end. The arguments are the same either
way, and the example below shows the paired-end one.



```
scout <- parameter_scout_paired_end(
  annotation_file          = "annotation.gff3",
  original_sRNA_annotation = "unknown",
  strandedness             = "reversely_stranded")   # or "stranded"
```



`"annotation.gff3"` is a placeholder. Replace it with your own GFF filename,
exactly as it appears in your folder.

**`strandedness` describes your library, not your preference.** Most bacterial
RNA-seq kits use a dUTP protocol, which is `"reversely_stranded"`. A forward
protocol is `"stranded"`, which is also the default if you leave the argument
out. Getting this wrong does not stop the scan and does not warn you: it reports
plausible numbers computed from the wrong strand. If you do not know which your
library is, the kit documentation will say, or RSeQC's `infer_experiment.py` will
tell you from a BAM file.

Give `original_sRNA_annotation` the same value you will later give `feature_file_editor()`, so both stages measure the same regions.

### Step 4\. Wait

The scan reads every BAM file once, which costs roughly the same as one `feature_file_editor()` pass. It is not quick on deep data and can take 3 minutes or more per BAM.

### Step 5\. Read the figure first, then the numbers

The figure is the main output and the numbers are secondary. It has two parts:

* a **table** of percentiles, one row per BAM file, with a bold **Median across BAMs** row at the bottom. The suggested cutoffs come from the bold row.  
* a **distribution** for each BAM file, drawn as a violin with a box plot below it, on a log2 scale, with dashed markers at P25, the median and P75 to make it easier to read.

On a shallow library one of the three markers may be missing. That happens when
the P25 cell in the bold row reads 1: `log2(1)` sits on the axis edge so no line
can be drawn there, and a cutoff of 1 would keep everything shown in any case. The
value is still in the table, and nothing is wrong.

Look at the shapes. Are the files similar to each other, or is one much deeper than the rest? Is there a long tail? Does the marker for the median sit on the busiest part of the distribution or on a quiet slope? The annotation block on the right gives you the two pairs, and the note beneath says what it should: try these first, other cutoffs may suit your data better.

### Step 6\. Record the coverage pairs, then type them in

Both suggested cutoff pairs are written to `scout_suggestions.txt`, in the `scout` folder alongside the other outputs.

**Record your chosen values**, then type them into `feature_file_editor()` by hand. Parameter Scout will not do this for you.

### The output files for Parameter Scout



| file | what it is |
| --- | --- |
| `scout_percentiles.tsv` | the percentile table, one row per BAM file |
| `scout_distribution.tsv` | the coverage values and their weights, split by file and strand |
| `scout_suggestions.txt` | both pairs, written out to be copied by hand |
| `scout_run_log.txt` | what was run: files, settings, R version, elapsed time |
| `scout_figure_<label>.png` | the figure |
| `scout_figure_<label>.pdf` | the same figure, for printing |



### If something goes wrong



| what you see | what it means |
| --- | --- |
| `argument "annotation_file" is missing, with no default`, or the same for `original_sRNA_annotation` | these are the only two you must always supply. Everything else has a default |
| `there is no package called ...` | that package is not installed. See the list above step 1 and install it |
| `could not find function "major_features"` | you skipped the second `source()` line in step 2 |
| `could not find function "parameter_scout_paired_end"`, or the same for `parameter_scout_single_end` | you loaded one of the two modules and called the other one's function. See step 3 |
| `Annotation file not found` | the GFF name is mistyped, or the file is not in your working directory |
| `No BAM files found` | wrong folder, or nothing in it ends in `.bam`. Upper or lower case both count |
| `Invalid 'strandedness' value` | you passed `"unstranded"`, which prediction cannot use |
| a message about scanning up to 9 BAM files | you gave it ten or more. Scan in subsets, naming them with `bam_txt_list` |
| `parameter scout expects a single reference sequence` | the BAM header lists more than one reference sequence. Parameter Scout handles single-chromosome genomes only |
| a warning that a file `yielded no reads after filtering`, and `NA` percentiles | one file means that file is empty or badly filtered. Every file means the filter is wrong: check `mapqFilter` against your aligner, and check you are using the module that matches your reads, since the paired-end one discards single-end reads entirely |
| the files are not where you expected | the working directory was not the folder you thought. The scan prints the full path of every file it writes, so check those and check `getwd()` |
| the numbers look wrong for the dataset | check the run log for which BAM files were actually read |



---

## References

Ozuna, A., Liberto, D., Joyce, R. M., Arnvig, K. B., & Nobeli, I. (2020). baerhunter: An R package for the discovery and analysis of expressed non-coding regions in bacterial RNA-seq data. _Bioinformatics_, _36_(3), 966-969. [https://doi.org/10.1093/bioinformatics/btz643](https://doi.org/10.1093/bioinformatics/btz643)

Sivasankaran, A. (2024). _Detection and annotation of non-coding RNAs in bacteria using baerhunter: Developing a user-friendly RShiny application with intergenic coverage analysis for informed parameter selection_ [Unpublished MSc dissertation]. Birkbeck, University of London.
