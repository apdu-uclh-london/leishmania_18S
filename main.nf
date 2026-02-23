nextflow.enable.dsl=2

// Parameters
params.in          = null
params.out         = "screening_results"
params.host_index  = null
params.ref_18s     = "${projectDir}/ref/leish_18S_ref.fasta"

// Includes
include { PREP_FASTQS }   from './modules/prep'
include { FASTP }         from './modules/fastp'
include { BOWTIE2_HRR }   from './modules/bowtie2'
include { ALIGN_18S }     from './modules/alignment_18S'
include { SCREEN_REPORT } from './modules/screen_report'

workflow {
    if (!params.in || !params.host_index) {
        error "Usage: nextflow run main.nf --in <dir> --host_index <path> --out <dir>"
    }

    // 1. Setup Inputs
    input_ch = Channel.fromFilePairs("${params.in}/*_R{1,2}*.fastq.gz")
    ref_ch = Channel.fromPath("${params.ref_18s}*").collect()
    host_ch = Channel.fromPath("${params.host_index}*").collect()

    // 2. Execution Flow
    PREP_FASTQS(input_ch)
    FASTP(PREP_FASTQS.out.merged_reads)
    BOWTIE2_HRR(FASTP.out.trimmed_reads, host_ch)
    ALIGN_18S(BOWTIE2_HRR.out.clean_reads, ref_ch)
    
    // 3. Reporting
    SCREEN_REPORT(ALIGN_18S.out.stats)
}