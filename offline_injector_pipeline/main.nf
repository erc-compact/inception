#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { filtool } from './processes'
include { fold_par } from './processes'
include { peasoup0 } from './processes'
include { peasoup1 } from './processes'
include { peasoup2 } from './processes'
include { peasoup3 } from './processes'
include { candidate_filter } from './processes'
include { fold_cand } from './processes'
include { score_collect } from './processes'
include { get_results } from './processes'



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
    inj_fold_par = fold_par(inj_pulsars)
    inj_peasoup = peasoup_spawner(inj_filtool)
    inj_cand_filter = candidate_filter(inj_peasoup)
    inj_fold_cand = fold_cand(inj_cand_filter)
    inj_results = score_collect(inj_fold_cand, inj_fold_par)

    emit:
        inj_results
}   

workflow {
    injection_batch = Channel.from(1..params.n_injections) 
    inj_results = injection_batch | injection_pipeline

    // get_results(0)
}

