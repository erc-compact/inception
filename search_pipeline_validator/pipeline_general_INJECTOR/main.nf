#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { pulsarx_parfold } from './processes'
include { presto_parfold } from './processes'
include { dspsr_parfold } from './processes'

workflow INJECT {
    take:
        injection_number

    main:
    inj_pulsars = injection(injection_number)
    inj_px = pulsarx_parfold(inj_pulsars)

    // inj_presto = presto_parfold(inj_pulsars)
    // inj_dspsr = dspsr_parfold(inj_pulsars)
}   

workflow {
    injection_batch = Channel.from(params.start..params.end) 
    inj_results = injection_batch | INJECT
}

