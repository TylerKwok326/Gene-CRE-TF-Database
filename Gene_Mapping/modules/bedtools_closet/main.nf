#!/usr/bin/env nextflow

process BEDTOOLS_CLOSEST {
    label 'process_single'
    publishDir "$params.outdir", mode: 'copy'
    container 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

    input: 
    path(tss)
    tuple val(name), path(cres)

    output: 
    tuple val(name), path("${name}.closest.bed")

    shell:
    """
    zcat $tss | sort -k1,1 -k2,2n > sorted_tss.bed
    tail -n +2 $cres | sort -k1,1 -k2,2n > sorted_cres.bed
    bedtools closest -a sorted_cres.bed -b sorted_tss.bed -d > ${name}.closest.bed
    """
}