nextflow.enable.dsl=2


process injection {
    maxForks params.batch_size
    label "injection"
    container params.injection_image

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
    container params.filtool_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}

    python3.6 ${params.pipeline_code}/pipeline_filtool.py  --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --threads=${task.cpus}

    """
}

process pulsarx_parfold {
    label "pulsarx_parfold"
    container params.fold_par_image

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


process pulsarx_candfold {
    label "pulsarx_candfold"
    container params.fold_cand_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}
    source ${params.dependencies_config} ${params.tmp_dir} python3.6

    python3.6 ${params.pipeline_code}/pipeline_pulsarx_candfold.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process candidate_filter {
    label "candidate_filter"
    container params.candidate_filter_image

    input:
        val injection_number

    output:
        val injection_number
    
    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}
    
    python3.6 ${params.pipeline_code}/MMGPS_candidate_filter.py --processing_args=${params.config_params} --out_dir=${params.output_dir} --injection_number=${injection_number}
    
    """
}

process pics_scorer {
    label "pics_scorer"
    container params.pics_scorer_image

    input:
        val inj_cand

    output:
        val inj_cand

    scratch params.tmp_dir

    script:
    """
    python3.6 ${params.pipeline_code}/MMGPS_PICS_scorer.py --processing_args=${params.processing_args} --injection_number=${inj_cand}  --out_dir=${params.output_dir}

    """
}
