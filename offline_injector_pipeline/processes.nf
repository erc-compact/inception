nextflow.enable.dsl=2

process injection_setup {
    label "injection_setup"

    input:
    val injection_number

    scratch "${params.injection_setup.tmp_dir}/inj_${injection_number}"

    script:
    """

    command_inj="python3 ${projectDir}/pipeline_injector.py"
    inputs_inj="--search_args=${params.search_params} --injection_file=${params.injection_plan} --data_dir=${params.data_dir} --out_dir=${params.output_dir} --injection_number=${injection_number} --ephem=${params.ephem} --ncpus=${task.cpus}"
    
    ${params.injection_setup.injection_image} $command_inj $inputs_inj


    command_fil="python3.6 ${projectDir}/pipeline_filtool.py"
    inputs_fil="--injection_number=${injection_number} --search_args=${params.search_params} --out_dir=${params.output_dir}  --ncpus=${task.cpus}"

    ${params.injection_setup.filtool_image} bash -c "
        export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
        export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
        export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
        $command_fil $inputs_fil"

    """
}

process fold_par {
    label "fold_par"

    input:
    val injection_number

    scratch "${params.fold_par.tmp_dir}/inj_${injection_number}"

    script:
    """

    command="python3.6 ${projectDir}/pipeline_fold.py"
    inputs="--mode=par --search_args=${params.search_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}"

    ${params.fold_par.image} $command $inputs
    """
}

process peasoup {
    label "peasoup"

    input:
    tuple val(injection_number), val(tscrunch)

    scratch "${params.peasoup.tmp_dir}/inj_${injection_number}_${tscrunch}"

    script:
    """

    command="python3.6 ${projectDir}/pipeline_peasoup.py"
    inputs="--tscrunch_index=${tscrunch} --search_args=${params.search_params} --out_dir=${params.output_dir}  --injection_number=${injection_number}"

    ${params.peasoup.image} bash -c "
        export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
        export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
        export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
        $command $inputs"

    """
}

process candidate_filter {
    label "candidate_filter"

    input:
    val injection_number

    scratch "${params.candidate_filter.tmp_dir}/inj_${injection_number}"

    script:
    """

    command="python3.6 ${projectDir}/pipeline_candidate_filter.py"
    inputs="--search_args=${params.search_params} --data_dir=${params.data_dir} --out_dir=${params.output_dir}  --injection_number=${injection_number}"

    ${params.candidate_filter.image} $command $inputs

    """
}

process fold_cand {
    label "fold_cand"

    input:
    val injection_number

    scratch "${params.fold_cand.tmp_dir}/inj_${injection_number}"

    script:
    """

    command="python3.6 ${projectDir}/pipeline_candidate_filter.py"
    inputs="--mode=cand --search_args=${params.search_params} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}"

    ${params.fold_cand.image} $command $inputs

    """
}