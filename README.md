# methseq 

- My Entire genome methylation pipeline wrote earlier for work at Poland. 
- An end to end methylation analysis. 
- Releasing now for public usage.

```
# Steps:
#   1. Raw read QC                (FastQC)
#   2. Adapter/quality trimming   (Trim Galore, bisulfite-aware)
#   3. Post-trim QC               (FastQC)
#   4. Genome preparation         (Bismark, run once, cached)
#   5. Bisulfite alignment        (Bismark / Bowtie2)
#   6. Deduplication              (deduplicate_bismark)
#   7. Methylation extraction     (bismark_methylation_extractor)
#   8. Genome-wide report         (bismark2report / bismark2summary)
#   9. Coverage -> CpG report     (coverage2cytosine)
#  10. Differential methylation   (methylKit, via Rscript)
```

Gaurav Sablok \
gsablok@proton.me
