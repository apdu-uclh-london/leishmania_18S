process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*.fq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*.json"), emit: json

    script:
    if (params.single_end) {
        """
        fastp \
            -i ${reads[0]} \
            -o ${sample_id}_trim.fq.gz \
            --length_required 50 \
            --thread ${task.cpus} \
            --json ${sample_id}.fastp.json
        """
    } else {
        """
        fastp \
            -i ${reads[0]} -I ${reads[1]} \
            -o ${sample_id}_trim_R1.fq.gz -O ${sample_id}_trim_R2.fq.gz \
            --detect_adapter_for_pe \
            --length_required 50 \
            --thread ${task.cpus} \
            --json ${sample_id}.fastp.json
        """
    }
}