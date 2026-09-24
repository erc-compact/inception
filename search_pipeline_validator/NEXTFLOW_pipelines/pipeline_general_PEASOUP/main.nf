#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { filtool } from './processes'
include { pulsarx_parfold } from './processes'
include { peasoup_search } from './processes'
include { peasoup_setup } from './processes'
include { match_candidates } from './processes'
include { pulsarx_candfold } from './processes'
include { classifier } from './processes'
include { collector } from './processes'



def expand_plan(channel) {
    channel.flatMap { injection_number, plan ->
        def segments = plan.readLines()
        def key = groupKey(injection_number, segments.size())

        segments.collect { segment ->
            tuple(key, segment.trim())
        }
    }
}

def collapse_tag(channel) {
    channel
        .groupTuple()
        .map { key, items ->
            key.getGroupTarget()
        }
}


workflow PEASOUP {
    take:
        cand_channel

    main:
        peasoup_jobs = expand_plan(peasoup_setup(cand_channel))

        search_jobs = collapse_tag(peasoup_search(peasoup_jobs))

    emit:
        search_jobs
}


workflow INJECT {
    take:
        injection_number

    main:
        inj_pulsars = injection(injection_number)

        inj_filtool = filtool(inj_pulsars)

        inj_fold_par = pulsarx_parfold(inj_pulsars)
        
        inj_peasoup = PEASOUP(inj_filtool)

        inj_tag_match = expand_plan(match_candidates(inj_peasoup))

        inj_fold_cand = pulsarx_candfold(inj_tag_match)

        inj_classified = classifier(inj_fold_cand)

        inj_output = collapse_tag(inj_classified)
    emit:
        inj_output

}

workflow {
    injection_batch = Channel.from(params.start..params.end) 

    inj_results = INJECT(injection_batch)

    collector(inj_results.collect())
}

