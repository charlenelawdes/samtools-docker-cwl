process MERGE {
    container="ghcr.io/bwbioinfo/samtools-docker-cwl:e80764711a121872e9ea35d90229cec6dd6d8dec"
    tag "merge $chunk"
    label "sam_mid"
    
    input:
    tuple val(sample_id), val(chunk), path(bam)
    val parameters
    
    output:
    tuple val(sample_id), path("${sample_id}_${chunk}.bam")

    script:
    """
    samtools merge \\
        -@ ${params.threads} \\
        $parameters \\
        "${sample_id}_${chunk}.bam" \\
        $bam
    """
}

process SORT {
    publishDir "${params.out_dir}/alignments", mode: 'link', enabled: params.publish
    label "sam_big"
    container="ghcr.io/bwbioinfo/samtools-docker-cwl:e80764711a121872e9ea35d90229cec6dd6d8dec"
    tag "sam_sort $sample_id"

    input: 
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sam.baseName}.bam"), path("${sam.baseName}.bam.bai")
    
    script:
    """
    samtools sort \\
    -@ $params.threads \\
    --write-index \\
    -o ${sam.baseName}.bam##idx##${sam.baseName}.bam.bai \
    $bam
    """
}