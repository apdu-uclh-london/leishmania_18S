process SCREEN_REPORT {
    publishDir "${params.out}/results", mode: 'copy'

    input:
    tuple val(sample_id), val(raw), val(trimmed), val(clean), val(mapped)

    output:
    path "${sample_id}_screen_results.csv", emit: csv

    script:
    """
    python3 <<EOF
    import csv
    
    count = int("${mapped}")
    
    if count > 500:
        result, status = "POSITIVE", "PASS"
    elif count > 50:
        result, status = "INCONCLUSIVE", "REVIEW"
    else:
        result, status = "NEGATIVE", "FAIL"

    with open("${sample_id}_screen_results.csv", 'w') as f:
        f.write("Sample_ID,Leishmania_Detected,Raw_Reads,Trimmed_Reads,Decontaminated_Reads,Total_Mapped_Reads,QC_Status\\n")
        f.write(f"${sample_id},{result},${raw},${trimmed},${clean},{count},{status}\\n")
    EOF
    """
}