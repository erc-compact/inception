nextflow.enable.dsl=2

process setup {
    label "setup"
    executor 'local'

    output:
        path("tags.txt")

    script:
    """

    python3 ${params.pipeline_code}/preprocess_setup.py  --processing_args=${params.pipeline_config} --filterbanks=${params.filterbanks} --out_dir=${params.output_dir}

    """
}


process filtool {
    label "filtool"
    container params.pulsarx_image

    input:
        val(tag)

    output:
        val(tag)

    scratch params.tmp_dir

    script:
    """

    python3 ${params.pipeline_code}/preprocess_filtool.py  --tag=${tag} --processing_args=${params.pipeline_config} --out_dir=${params.output_dir} --threads=${task.cpus}

    """
}

process nullsar_zap {
    label "nullsar_zap"
    container params.nullsar_image

    input:
        val(tag) 
        val(mode)
        
    output:
        val(tag)

    scratch params.tmp_dir

    script:
    """

    python3 ${params.pipeline_code}/nullsar_zapper.py  --tag=${tag} --mode=${mode} --processing_args=${params.pipeline_config} --out_dir=${params.output_dir} --ncpus=${task.cpus} 

    """
}

process nullsar_fold {
    label "nullsar_fold"
    container params.pulsarx_image

    input:
        val(tag) 
        val(mode)

    output:
        val(tag)

    scratch params.tmp_dir

    script:
    """

    python3 ${params.pipeline_code}/nullsar_folder.py --tag=${tag} --mode=${mode} --processing_args=${params.pipeline_config} --out_dir=${params.output_dir} --ncpus=${task.cpus} 

    """
}


