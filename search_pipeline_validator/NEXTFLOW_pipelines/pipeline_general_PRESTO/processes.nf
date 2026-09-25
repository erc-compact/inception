nextflow.enable.dsl=2


process injection {
    maxForks params.batch_size
    label "injection"
    container params.python_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_general/pipeline_injector.py --processing_args=${params.config_params} --injection_plan=${params.injection_plan} --out_dir=${params.output_dir} --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process presto_parfold {
    label "presto_parfold"
    container params.presto_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_parfold.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process rfifind {
    label "rfifind"
    container params.presto_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_rfifind.py  --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --threads=${task.cpus}

    """
}

process presto_ddplan_setup {
    label "presto_setup"
    container params.python_image

    input:
        val injection_number

    output:
        tuple val(injection_number), path("*_DDPLAN_PLAN.txt")

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_setup.py --mode=ddplan --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process presto_search_setup {
    label "presto_setup"
    container params.python_image

    input:
        val injection_number

    output:
        tuple val(injection_number), path("*_PROCESS_PLAN.txt")

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_setup.py --mode=search --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process presto_dedisperse {
    label "presto_dedisperse"
    container params.presto_image

    input:
        tuple val(injection_number), val(ddplan)

    output:
        tuple val(injection_number), val(ddplan)

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_dedisperse.py --tag=${ddplan} --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process presto_fft {
    label "presto_fft"
    container params.presto_image

    input:
        tuple val(injection_number), val(segment)

    output:
        tuple val(injection_number), val(segment)

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_fft.py --tag=${segment} --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process presto_accelsearch {
    label "presto_accelsearch"
    container params.presto_image

    input:
        tuple val(injection_number), val(segment)

    output:
        tuple val(injection_number), val(segment)

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_accelsearch.py --tag=${segment} --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}


process presto_cleanup {
    label "presto_cleanup"
    container params.python_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_cleanup.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process presto_sift {
    label "presto_sift"
    container params.presto_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_scripts/ACCEL_sift.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process match_candidates {
    label "match_candidates"
    container params.python_image

    input:
        val injection_number

    output:
        tuple val(injection_number), path("*_FOLD_PLAN.txt")

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_general/pipeline_candidate_matcher.py --processing_args=${params.config_params} --out_dir=${params.output_dir} --injection_number=${injection_number}

    """
}

process presto_candfold {
    label "presto_candfold"
    container params.presto_image

    input:
        tuple val(injection_number), val(segment)

    output:
        tuple val(injection_number), val(segment)

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_presto/pipeline_presto_candfold.py --tag=${segment} --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process collector {
    label "collector"
    container params.python_image

    input:
        val results

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/PY_general/pipeline_collect_results.py --processing_args=${params.config_params} --out_dir=${params.output_dir}

    """
}
