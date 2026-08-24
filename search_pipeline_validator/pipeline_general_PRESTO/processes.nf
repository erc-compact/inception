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
    python3 ${params.pipeline_code}/pipeline_presto_parfold.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

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
    python3 ${params.pipeline_code}/pipeline_presto_rfifind.py  --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --threads=${task.cpus}

    """
}

process presto_search {
    label "presto_search"
    container params.presto_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    python3 ${params.pipeline_code}/pipeline_presto_search.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --threads=${task.cpus}

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
    python3 ${params.pipeline_code}/scripts/ACCEL_sift.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process presto_candfold {
    label "presto_candfold"
    container params.presto_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir
    
    script:
    """
    python3 ${params.pipeline_code}/pipeline_presto_candfold.py --processing_args=${params.config_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --threads=${task.cpus}

    """
}
