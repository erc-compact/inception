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
    python3 ${projectDir}/pipeline_injector.py --search_args=${params.search_params} --injection_file=${params.injection_plan} --data_dir=${params.data_dir} --out_dir=${params.output_dir} --injection_number=${injection_number} --ephem=${params.ephem} --ncpus=${task.cpus}
    
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

    python3.6 ${projectDir}/pipeline_filtool.py --injection_number=${injection_number} --search_args=${params.search_params} --out_dir=${params.output_dir}  --ncpus=${task.cpus}

    """
}

process fold_par {
    label "fold_par"
    container params.fold_par_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}
    source ${params.dependencies_config} ${params.tmp_dir}

    python3.6 ${projectDir}/pipeline_fold.py --mode=par --search_args=${params.search_params}  --injection_file=${params.injection_plan} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}


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

    python3.6 ${projectDir}/pipeline_peasoup.py --tscrunch_index=0 --search_args=${params.search_params} --injection_file=${params.injection_plan} --out_dir=${params.output_dir} --data_dir=${params.data_dir} --injection_number=${injection_number}

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

    python3.6 ${projectDir}/pipeline_peasoup.py --tscrunch_index=1 --search_args=${params.search_params} --injection_file=${params.injection_plan} --out_dir=${params.output_dir} --data_dir=${params.data_dir} --injection_number=${injection_number}

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

    python3.6 ${projectDir}/pipeline_peasoup.py --tscrunch_index=2 --search_args=${params.search_params} --injection_file=${params.injection_plan} --out_dir=${params.output_dir} --data_dir=${params.data_dir} --injection_number=${injection_number}

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

    python3.6 ${projectDir}/pipeline_peasoup.py --tscrunch_index=3 --search_args=${params.search_params} --injection_file=${params.injection_plan} --out_dir=${params.output_dir} --data_dir=${params.data_dir} --injection_number=${injection_number}

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
    
    python3.6 ${projectDir}/pipeline_candidate_filter.py --search_args=${params.search_params} --data_dir=${params.data_dir} --out_dir=${params.output_dir} --injection_number=${injection_number}
    
    """
}

process fold_cand {
    label "fold_cand"
    container params.fold_cand_image

    input:
        val injection_number

    output:
        val injection_number

    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}
    source ${params.dependencies_config} ${params.tmp_dir}

    python3.6 ${projectDir}/pipeline_fold.py --mode=cand --search_args=${params.search_params}  --injection_file=${params.injection_plan} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process pics_scorer {
    label "pics_scorer"
    container params.pics_scorer_image

    input:
        val inj_cand
        val inj_fold

    output:
        val inj_cand

    scratch params.tmp_dir

    script:
    """
    source ${params.singularity_config}

    python3.6 ${projectDir}/pipeline_PICS.py --pics_code=${params.pics_code} --pics_models=${params.pics_models} --injection_number=${inj_cand}  --out_dir=${params.output_dir}

    """
}