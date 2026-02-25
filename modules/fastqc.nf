process FASTQC {
    tag "$sample_id"
    publishDir "${params.out}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
    """
}