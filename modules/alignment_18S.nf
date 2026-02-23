process ALIGN_18S {
    tag "$sample_id"
    publishDir "${params.out}/alignment_18s", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path ref_18s

    output:
    tuple val(sample_id), env(TOTAL_MAPPED), emit: stats
    tuple val(sample_id), path("${sample_id}_18s.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} ${ref_18s} ${reads[0]} ${reads[1]} | \\
    samtools view -b | samtools sort -o ${sample_id}_18s.bam
    
    # Get count of mapped reads
    COUNT=\$(samtools view -c -F 4 ${sample_id}_18s.bam)
    echo "TOTAL_MAPPED=\$COUNT" > count.env
    source count.env
    """
}