#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { setup } from './processes'
include { filtool } from './processes'
include { nullsar_fold } from './processes'
include { nullsar_zap } from './processes'



workflow nullsar {
    take:
        tag_file

    main:

    fold_out_0 = nullsar_fold(tag_file, "INIT")
    zap_out_0 = nullsar_zap(fold_out_0, "INIT")

    fold_out_1 = nullsar_fold(zap_out_0, "OPTIMISE")
    zap_out_1 = nullsar_zap(fold_out_1, "NULL")

    fold_out_2 = nullsar_fold(zap_out_1, "CONFIRM")

    emit:
       zap_out_1

}


workflow {
    tag_file = setup()

    beam_channel = tag_file
        .splitText()
        .map { line -> line.trim() }

    filtool_out = filtool(beam_channel)

    nullsar_out = nullsar(filtool_out)

}