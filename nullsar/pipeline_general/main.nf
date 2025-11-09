#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { setup } from './processes'
include { filtool } from './processes'

include { nullsar_fold as nullsar_fold_init } from './processes'
include { nullsar_fold as nullsar_fold_opt }  from './processes'
include { nullsar_fold as nullsar_fold_conf } from './processes'

include { nullsar_zap  as nullsar_zap_init }  from './processes'
include { nullsar_zap  as nullsar_zap_null }  from './processes'



workflow nullsar {
    take:
        beam_channel

    main:

    fold_out_0 = nullsar_fold_init(beam_channel, "INIT")
    zap_out_0 = nullsar_zap_init(fold_out_0, "INIT")

    fold_out_1 = nullsar_fold_opt(zap_out_0, "OPTIMISE")
    zap_out_1 = nullsar_zap_null(fold_out_1, "NULL")

    fold_out_2 = nullsar_fold_conf(zap_out_1, "CONFIRM")

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