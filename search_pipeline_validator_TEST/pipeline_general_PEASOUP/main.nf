#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { filtool } from './processes'
include { pulsarx_parfold } from './processes'
include { peasoup_search } from './processes'
include { peasoup_setup } from './processes'
include { match_candidates } from './processes'
include { pulsarx_candfold } from './processes'


workflow PEASOUP {
    take: 
        cand_channel

    main:
        peasoup_jobs = peasoup_setup(cand_channel)
            .flatMap { tag, plan ->

                plan.readLines()
                    .collect { seg -> tuple(tag, seg.trim()) }
            }

        search_jobs = peasoup_search(peasoup_jobs)
        collapsed = search_jobs
            .map { tag, result -> tag }
            .distinct()

    emit:
        collapsed
}


workflow MATCHER {
    take: 
        cand_channel

    main:

        inj_match = match_candidates(cand_channel)
                .flatMap { tag, plan ->

                plan.readLines()
                    .collect { seg -> tuple(tag, seg.trim()) }
            }

    emit:
        inj_match

}


workflow injection_pipeline {
    take:
        injection_number

    main:
    inj_pulsars = injection(injection_number)

    inj_filtool = filtool(inj_pulsars)
    inj_fold_par = pulsarx_parfold(inj_pulsars)

    inj_peasoup = PEASOUP(inj_filtool)

    inj_match = MATCHER(inj_peasoup)

    inj_fold_cand = pulsarx_candfold(inj_match)
}

workflow {
    injection_batch = Channel.from(params.start..params.end) 
    inj_results = injection_batch | injection_pipeline
}

