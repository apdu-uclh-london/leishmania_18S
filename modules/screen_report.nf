process SCREEN_REPORT {
    publishDir "${params.out}/results", mode: 'copy'

    input:
    tuple val(sample_id), val(mapped_count)

    output:
    path "${sample_id}_screen_results.csv"

    script:
    """
    python3 <<EOF
    import csv
    
    count = int("${mapped_count}")
    
    if count > 500:
        result = "POSITIVE"
        status = "PASS"
    elif count > 50:
        result = "INCONCLUSIVE"
        status = "REVIEW"
    else:
        result = "NEGATIVE"
        status = "PASS"

    with open("${sample_id}_screen_results.csv", 'w') as f:
        f.write("Sample_ID,Leishmania_Detected,Mapped_18S_Reads,QC_Status\\n")
        f.write(f"${sample_id},{result},{count},{status}\\n")
    EOF
    """
}