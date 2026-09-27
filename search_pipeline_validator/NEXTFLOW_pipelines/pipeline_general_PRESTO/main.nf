#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { injection } from './processes'
include { presto_parfold } from './processes'
include { rfifind } from './processes'
include { filtool } from './processes'
include { presto_ddplan_setup } from './processes'
include { presto_search_setup } from './processes'
include { presto_dedisperse } from './processes'
include { presto_fft } from './processes'
include { presto_accelsearch } from './processes'
include { presto_cleanup } from './processes'
include { presto_sift } from './processes'
include { match_candidates } from './processes'
include { presto_candfold } from './processes'
include { classifier } from './processes'
include { collector } from './processes'



def presto_rfi_cleaner() {
    def config = new groovy.json.JsonSlurper().parseText(file(params.config_params).text)
    def search = config.presto_search_args ?: [:]
    def cleaner = search.rfi_cleaner ?: 'rfifind'

    if (!(cleaner in ['rfifind', 'filtool'])) {
        error "presto_search_args.rfi_cleaner must be 'rfifind' or 'filtool', not '${cleaner}'"
    }

    if (cleaner == 'filtool') {
        def conflicts = []
        if (search.mask) conflicts << 'presto_search_args.mask'
        if (config.presto_candfold_args?.mask) conflicts << 'presto_candfold_args.mask'
        if (search.birdies == 'rfifind') conflicts << 'presto_search_args.birdies'
        if (config.presto_parfold_args?.mask == 'rfifind') conflicts << 'presto_parfold_args.mask'
        if (conflicts) {
            error "rfi_cleaner is 'filtool', so ${conflicts.join(', ')} must not use rfifind or a mask - only one RFI cleaner can be used"
        }

        def f_args = config.filtool_args
        if (!f_args) error "rfi_cleaner is 'filtool' but the config has no filtool_args"
        if (!(1 in f_args.tscrunch)) error "filtool_args.tscrunch must include 1: PRESTO dedisperses the full-resolution filtool output"
        if (!f_args.save_filtool_fb) error "rfi_cleaner is 'filtool', so filtool_args.save_filtool_fb must be true"
    }

    return cleaner
}


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
        cleaned_channel

    main:
        ddplan_jobs = expand_plan(presto_ddplan_setup(cleaned_channel))
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

        if (presto_rfi_cleaner() == 'filtool') {
            inj_cleaned = filtool(inj_pulsars)
        } else {
            inj_cleaned = rfifind(inj_pulsars)
        }

        inj_fold_par = presto_parfold(inj_pulsars)

        inj_search = PRESTO(inj_cleaned)

        inj_sift = presto_sift(inj_search)

        inj_tag_match = expand_plan(match_candidates(inj_sift))

        inj_fold_cand = presto_candfold(inj_tag_match)

        inj_classified = classifier(inj_fold_cand)

        inj_output = collapse_tag(inj_classified)
    emit:
        inj_output

}

workflow {
    injection_batch = Channel.from(params.start..params.end)

    inj_results = INJECT(injection_batch)

    collector(inj_results.toList())
}
