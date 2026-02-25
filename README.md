# Leishmania 18S Amplicon Screening Pipeline

Nextflow DSL2 pipeline for clinical screening of Leishmania 18S amplicon sequencing data.

## Prerequisites
* Nextflow (v22.10.0+)
* Docker (uses `mbyott/leishmania_speciation:v1.2` container)
* Bowtie2 host index (e.g., T2T-CHM13v2)
* BWA 18S reference index (`ref/18s_ref.fasta`)

## Execution
Run the following command from the project root:

```bash
nextflow run main.nf \
  --in "/path/to/raw/FASTQ/dir" \
  --out "/path/to/output/dir" \
  --host_index "/path/to/db/T2T-CHM13v2" \
  --single_end false \
  -profile docker
```

Diagnostic Thresholds

Evaluated on ```Total_Mapped_Reads```:

*  500+: POSITIVE (QC: PASS)
*  51 - 500: INCONCLUSIVE (QC: REVIEW)
*  <= 50: NEGATIVE (QC: FAIL)
