#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { filtool } from './processes'
include { pulsarx_parfold } from './processes'
include { pulsarx_candfold } from './processes'
include { candidate_filter } from './processes'
include { pics_scorer } from './processes'

include { peasoup0 } from './peasoup_processes'
include { peasoup1 } from './peasoup_processes'
include { peasoup2 } from './peasoup_processes'
include { peasoup3 } from './peasoup_processes'


workflow peasoup_spawner {
    take:
        injection_number

    main:

    p0 = peasoup0(injection_number)
    p1 = peasoup1(injection_number)
    p2 = peasoup2(injection_number)
    p3 = peasoup3(injection_number)

    emit:
       p0

}


workflow injection_pipeline {
    take:
        injection_number

    main:
    inj_pulsars = injection(injection_number)
    inj_filtool = filtool(inj_pulsars)
    inj_fold_par = pulsarx_parfold(inj_pulsars)
    inj_peasoup = peasoup_spawner(inj_filtool)
    inj_filter = candidate_filter(inj_peasoup)
    inj_fold_cand = pulsarx_candfold(inj_peasoup)
    inj_pics = pics_scorer(inj_fold_cand)
}

workflow {
    injection_batch = Channel.from(1..params.n_injections) 
    inj_results = injection_batch | injection_pipeline
}

