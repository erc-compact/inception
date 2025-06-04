#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { pulsarx_parfold } from './processes'



workflow injection_pipeline {
    take:
        injection_number

    main:
    inj_pulsars = injection(injection_number)
    inj_fold_par = pulsarx_parfold(inj_pulsars)
}   

workflow {
    injection_batch = Channel.from(1..params.n_injections) 
    inj_results = injection_batch | injection_pipeline
}

