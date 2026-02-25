nextflow.enable.dsl=2

// Parameters
params.in          = null
params.out         = "screening_results"
params.host_index  = null
params.ref_18s     = "${projectDir}/ref/18s_ref.fasta"
params.single_end  = false

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

    // 1. Setup Inputs - Safely strip all Illumina extensions & lane IDs
    input_ch = Channel.fromPath("${params.in}/{,**/}*.{fastq,fq}.gz")
        .map { file -> 
            def id = file.name
                .replaceAll(/\.fastq\.gz$|\.fq\.gz$/, "") 
                .replaceAll(/_001$/, "")
                .replaceAll(/_R[12]$/, "")
                .replaceAll(/_[12]$/, "")
                .replaceAll(/_L\d{3}$/, "") 
                .replaceAll(/_S\d+$/, "")
            return tuple(id, file) 
        }
        .groupTuple()

    ref_ch = Channel.fromPath("${params.ref_18s}*").collect()
    host_ch = Channel.fromPath("${params.host_index}*").collect()

    // 2. Execution Flow
    PREP_FASTQS(input_ch)
    FASTP(PREP_FASTQS.out.merged_reads)
    BOWTIE2_HRR(FASTP.out.trimmed_reads, host_ch)
    ALIGN_18S(BOWTIE2_HRR.out.clean_reads, ref_ch)
    
    // 3. Reporting
    report_ch = PREP_FASTQS.out.raw_count
        .join(FASTP.out.trimmed_count)
        .join(BOWTIE2_HRR.out.read_count)
        .join(ALIGN_18S.out.stats)

    SCREEN_REPORT(report_ch)

    SCREEN_REPORT.out.csv
        .collectFile(
            name: 'master_18S_screen_results.csv', 
            keepHeader: true, 
            storeDir: "${params.out}/results"
        )
}