nextflow.enable.dsl=2


process peasoup0 {
    label "peasoup"
    container params.peasoup_image

    input:
        val injection_number
        
    output:
        val injection_number
        
    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}
    source ${params.dependencies_config} ${params.tmp_dir}

    python3.6 ${params.pipeline_code}/pipeline_peasoup.py --tscrunch_index=0 --processing_args=${params.search_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process peasoup1 {
    label "peasoup"
    container params.peasoup_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}
    source ${params.dependencies_config} ${params.tmp_dir}

    python3.6 ${params.pipeline_code}/pipeline_peasoup.py --tscrunch_index=1 --processing_args=${params.search_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process peasoup2 {
    label "peasoup"
    container params.peasoup_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}
    source ${params.dependencies_config} ${params.tmp_dir}

    python3.6 ${params.pipeline_code}/pipeline_peasoup.py --tscrunch_index=2 --processing_args=${params.search_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}

process peasoup3 {
    label "peasoup"
    container params.peasoup_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}
    source ${params.dependencies_config} ${params.tmp_dir}

    python3.6 ${params.pipeline_code}/pipeline_peasoup.py --tscrunch_index=3 --processing_args=${params.search_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}

    """
}