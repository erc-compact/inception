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
    python3 ${params.pipeline_code}/pipeline_injector.py --processing_args=${params.config_params} --injection_plan=${params.injection_plan} --out_dir=${params.output_dir} --injection_number=${injection_number} --ncpus=${task.cpus}
    
    """
}

process filtool {
    label "filtool"
    container params.pulsarx_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """

    python3 ${params.pipeline_code}/pipeline_filtool.py  --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --threads=${task.cpus}

    """
}

process pulsarx_parfold {
    label "pulsarx_parfold"
    container params.pulsarx_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """

    python3 ${params.pipeline_code}/pipeline_pulsarx_parfold.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process peasoup_setup {
    label "peasoup_setup"
    container params.python_image

    input:
        val injection_number
        
    output:
        tuple val(injection_number), path("*_PROCESS_PLAN.txt")
        
    scratch params.tmp_dir

    script:
    """

    python3 ${params.pipeline_code}/pipeline_peasoup_setup.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process peasoup_search {
    label "peasoup_search"
    container params.peasoup_image

    input:
        tuple val(injection_number), val(segment)
        
    output:
        tuple val(injection_number), val(segment)
        
    scratch params.tmp_dir

    script:
    """

    python3 ${params.pipeline_code}/pipeline_peasoup_search.py --tag=${segment} --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process pulsarx_candfold {
    label "pulsarx_candfold"
    container params.pulsarx_image

    input:
        tuple val(injection_number), val(segment)

    output:
        tuple val(injection_number), val(segment)

    scratch params.tmp_dir

    script:
    """

    python3 ${params.pipeline_code}/pipeline_pulsarx_candfold.py --tag=${segment} --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

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
    python3 ${params.pipeline_code}/pipeline_candidate_matcher.py --processing_args=${params.config_params} --out_dir=${params.output_dir} --injection_number=${injection_number}
    
    """
}