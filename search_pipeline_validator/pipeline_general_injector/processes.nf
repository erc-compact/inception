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
    source ${params.dependencies_config} ${params.tmp_dir}
    python3 ${params.pipeline_code}/pipeline_pulsarx_parfold.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

