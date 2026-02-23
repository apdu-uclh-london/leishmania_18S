process BOWTIE2_HRR {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/decontaminated", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}_final*.fq.gz"), emit: clean_reads
    path "${sample_id}_host_stats.txt", emit: stats
    tuple val(sample_id), env(CLEAN_COUNT), emit: read_count

    script:
    def idx_base = index_files.find { it.name.endsWith('.1.bt2') }.name.minus('.1.bt2')
    
    if (params.single_end) {
        """
        bowtie2 --end-to-end --sensitive -p ${task.cpus} -x ${idx_base} \
            -U ${reads[0]} -S output.sam 2> ${sample_id}.log

        samtools fastq -f 4 output.sam | gzip > ${sample_id}_final.fq.gz
                
        LINES=\$(zcat ${sample_id}_final.fq.gz | wc -l)
        CNT=\$((LINES / 4))
        
        echo "CLEAN_COUNT=\$CNT" > count.env
        source count.env
        cp ${sample_id}.log ${sample_id}_host_stats.txt
        """
    } else {
        """
        bowtie2 --end-to-end --sensitive -p ${task.cpus} -x ${idx_base} \
            -1 ${reads[0]} -2 ${reads[1]} -S output.sam 2> ${sample_id}.log

        samtools fastq -f 12 -1 ${sample_id}_final_1.fq.gz -2 ${sample_id}_final_2.fq.gz output.sam
        
        LINES=\$(zcat ${sample_id}_final_1.fq.gz | wc -l)
        CNT=\$((LINES / 2))
        
        echo "CLEAN_COUNT=\$CNT" > count.env
        source count.env
        cp ${sample_id}.log ${sample_id}_host_stats.txt
        """
    }
}