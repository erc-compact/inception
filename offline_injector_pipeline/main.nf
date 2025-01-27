#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { filtool } from './processes'
include { fold_par } from './processes'
include { peasoup } from './processes'
include { candidate_filter } from './processes'
include { fold_cand } from './processes'
include { score_collect } from './processes'
include { get_results } from './processes'


workflow injection_pipeline {
    take:
        injection_number

    main:
    inj_pulsars = injection(injection_number)
    inj_filtool = filtool(inj_pulsars)
    inj_fold_par = fold_par(inj_pulsars)
    // inj_filtool.view()
    // Split channels based on the no. of splits
    split_inj_filtool = inj_filtool
        .flatMap { injection_number, filterbanks -> 
            // Make sure that filterbanks is a list
            filterbank_list = filterbanks instanceof List ? filterbanks : [filterbanks]

            // Count the no of filterbanks
            n_filterbanks = filterbank_list.size()
        
            // Collect tuples in a list
            def tuples = []
            for (index in 0..<n_filterbanks) {
                tuples << tuple(injection_number, index, filterbank_list)
            }

        return tuples
    }
    
    // split_inj_filtool.view()
    inj_peasoup = peasoup(split_inj_filtool)
    inj_cand_filter = candidate_filter(inj_peasoup)
    inj_fold_cand = fold_cand(inj_cand_filter)
    inj_results = score_collect(inj_fold_cand, inj_fold_par)

    // emit:
    //     inj_results
}   

workflow {
    injection_batch = Channel.from(1..params.n_injections) 
    inj_results = injection_batch | injection_pipeline

    // get_results(0)
}

