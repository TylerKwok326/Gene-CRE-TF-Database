#!/usr/bin/env nextflow

include { BEDTOOLS_CLOSEST } from './modules/bedtools_closest'

workflow {

    Channel.fromPath(params.cres)
        .splitCsv(header: true)
        .map { row -> tuple(row.name, file(row.path)) }
        .set { cres_ch }

    BEDTOOLS_CLOSEST(params.tss, cres_ch)

}