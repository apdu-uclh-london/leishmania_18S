process PREP_FASTQS {
    tag "$sample_id"
    label 'process_low'
    
    // 1. SAVE OUTPUT: 'results/raw_fastq'
    publishDir "${params.outdir}/raw_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    // Captures "_merged.fastq.gz" (SE) OR "_R1_merged.../_R2_merged..." (PE)
    tuple val(sample_id), path("${sample_id}*_merged.fastq.gz"), emit: merged_reads

    script:
    if (params.single_end) {
        """
        echo "SE MODE: Merging lanes for ${sample_id}"
        
        # 1. List files, Sort by name (preserves Lane order), Merge to temp file
        # (Writing to .tmp first prevents 'input file is output file' errors)
        ls *fastq.gz | sort | xargs cat > merged.tmp
        
        # 2. Rename to final SE output
        mv merged.tmp ${sample_id}_merged.fastq.gz
        """
    } else {
        """
        echo "PE MODE: Merging lanes for ${sample_id}"
        
        # 1. Merge R1 files
        # We look for files containing 'R1' in the name
        ls *R1* | sort | xargs cat > r1.tmp
        mv r1.tmp ${sample_id}_R1_merged.fastq.gz
        
        # 2. Merge R2 files
        ls *R2* | sort | xargs cat > r2.tmp
        mv r2.tmp ${sample_id}_R2_merged.fastq.gz
        """
    }
}