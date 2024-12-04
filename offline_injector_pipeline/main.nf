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


workflow injection_pipeline {
    take:
    injection_number

    main:
    c0 = injection(injection_number)
    c1 = filtool(c0)
    c2 = fold_par(c1)
    c3 = peasoup3(c2)
    c4 = peasoup2(c3)
    c5 = peasoup1(c4)
    c6 = peasoup0(c5)
    c7 = candidate_filter(c6)
    c8 = fold_cand(c7)
}

workflow {

    injection_range = Channel.from(1..params.n_injections)

    injection_pipeline(injection_range)
    
}