process PREP_FASTQS {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.out}/raw_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*_merged.fastq.gz"), emit: merged_reads
    tuple val(sample_id), env(RAW_COUNT), emit: raw_count

    script:
    if (params.single_end) {
        """
        cat \$(ls *.fastq.gz | sort) > ${sample_id}_merged.fastq.gz
        LINES=\$(zcat ${sample_id}_merged.fastq.gz | wc -l)
        export RAW_COUNT=\$((LINES / 4))
        """
    } else {
        """
        cat \$(ls *R1*.fastq.gz | sort) > ${sample_id}_R1_merged.fastq.gz
        cat \$(ls *R2*.fastq.gz | sort) > ${sample_id}_R2_merged.fastq.gz
        LINES=\$(zcat ${sample_id}_R1_merged.fastq.gz | wc -l)
        export RAW_COUNT=\$((LINES / 2))
        """
    }
}