#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

    def instances = Channel.from(1..params.n_injections) 

    instances
        .map { injection_number -> workflow_instance(injection_number) }
}

workflow_instance(injection_number) {

    injection_setup(injection_number)

    fold_par(injection_number) &

    def peasoup_channel = Channel
        .from(0..(params.n_downsamp - 1)) 
        .map { tscrunch -> [injection_number, tscrunch] } 

    peasoup(peasoup_channel)

    candidate_filter(injection_number)

    fold_cand(injection_number)
}