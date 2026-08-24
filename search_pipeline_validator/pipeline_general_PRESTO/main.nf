#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { presto_parfold } from './processes'
include { rfifind } from './processes'
include { presto_search } from './processes'
include { presto_sift } from './processes'
include { presto_fold } from './processes'


workflow INJECT {
    take:
        injection_number

    main:
    inj_pulsars = injection(injection_number)

    inj_rfifind = rfifind(inj_pulsars)

    inj_fold_par = presto_parfold(inj_pulsars)

    inj_presto = presto_search(inj_rfifind)

    inj_sift = presto_sift(inj_presto)
    
    inj_fold = presto_fold(inj_sift)
}

workflow {
    injection_batch = Channel.from(params.start..params.end) 
    inj_results = injection_batch | INJECT
}
