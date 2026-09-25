#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { presto_parfold } from './processes'
include { rfifind } from './processes'
include { presto_ddplan_setup } from './processes'
include { presto_search_setup } from './processes'
include { presto_dedisperse } from './processes'
include { presto_fft } from './processes'
include { presto_accelsearch } from './processes'
include { presto_cleanup } from './processes'
include { presto_sift } from './processes'
include { match_candidates } from './processes'
include { presto_candfold } from './processes'
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


workflow PRESTO {
    take:
        rfifind_channel

    main:
        ddplan_jobs = expand_plan(presto_ddplan_setup(rfifind_channel))
        dedisp_done = collapse_tag(presto_dedisperse(ddplan_jobs))

        segment_jobs = expand_plan(presto_search_setup(dedisp_done))

        fft_jobs = presto_fft(segment_jobs)
        search_jobs = presto_cleanup(collapse_tag(presto_accelsearch(fft_jobs)))

    emit:
        search_jobs
}


workflow INJECT {
    take:
        injection_number

    main:
        inj_pulsars = injection(injection_number)

        inj_rfifind = rfifind(inj_pulsars)

        inj_fold_par = presto_parfold(inj_pulsars)

        inj_search = PRESTO(inj_rfifind)

        inj_sift = presto_sift(inj_search)

        inj_tag_match = expand_plan(match_candidates(inj_sift))

        inj_fold_cand = presto_candfold(inj_tag_match)

        inj_output = collapse_tag(inj_fold_cand)
    emit:
        inj_output

}

workflow {
    injection_batch = Channel.from(params.start..params.end)

    inj_results = INJECT(injection_batch)

    collector(inj_results.collect())
}
