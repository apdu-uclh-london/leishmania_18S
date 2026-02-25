process FASTP {
    tag "$sample_id"
    publishDir "${params.out}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trim_R*.fq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*.json"), emit: json
    tuple val(sample_id), env(TRIMMED_COUNT), emit: trimmed_count

    script:
    """
    fastp \\
        -i ${reads[0]} ${params.single_end ? '' : "-I ${reads[1]}"} \\
        -o ${sample_id}_trim_R1.fq.gz ${params.single_end ? '' : "-O ${sample_id}_trim_R2.fq.gz"} \\
        --detect_adapter_for_pe --length_required 50 --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json

    COUNT=\$(grep -A 10 '"after_filtering"' ${sample_id}.fastp.json | grep -m 1 '"total_reads"' | sed 's/[^0-9]//g')
    export TRIMMED_COUNT=\$COUNT
    """
}