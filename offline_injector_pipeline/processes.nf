nextflow.enable.dsl=2

process injection {
    label "injection"
    container "${params.injection_image}"

    input:
    val injection_number

    output:
    val injection_number

    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel

    python3 ${projectDir}/pipeline_injector.py --search_args=${params.search_params} --injection_file=${params.injection_plan} --data_dir=${params.data_dir} --out_dir=${params.output_dir} --injection_number=${injection_number} --ephem=${params.ephem} --ncpus=${task.cpus}
    
    """
}

process filtool {
    label "filtool"
    container "${params.filtool_image}"

    input:
    val injection_number

    output:
    val injection_number

    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel

    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
    python3.6 ${projectDir}/pipeline_filtool.py --injection_number=${injection_number} --search_args=${params.search_params} --out_dir=${params.output_dir}  --ncpus=${task.cpus}

    """
}

process fold_par {
    label "fold_par"
    container "${params.fold_par_image}"

    input:
    val injection_number

    output:
    val injection_number

    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;

    python3.6 ${projectDir}/pipeline_fold.py --mode=par --search_args=${params.search_params}  --injection_file=${params.injection_plan} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}

process peasoup0 {
    label "peasoup0"
    container "${params.peasoup_image}"

    input:
    val injection_number
    // tuple val(injection_number), val(tscrunch)
    // val tscrunch

    output:
    val injection_number

    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel

    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
    python3.6 ${projectDir}/pipeline_peasoup.py --tscrunch_index=0 --search_args=${params.search_params} --injection_file=${params.injection_plan} --out_dir=${params.output_dir}  --injection_number=${injection_number}
    inputs=""

    """
}

process peasoup1 {
    label "peasoup1"
    container "${params.peasoup_image}"

    input:
    val injection_number
    // tuple val(injection_number), val(tscrunch)
    // val tscrunch

    output:
    val injection_number

    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel

    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
    python3.6 ${projectDir}/pipeline_peasoup.py --tscrunch_index=1 --search_args=${params.search_params} --injection_file=${params.injection_plan} --out_dir=${params.output_dir}  --injection_number=${injection_number}
    inputs=""

    """
}

process peasoup2 {
    label "peasoup2"
    container "${params.peasoup_image}"

    input:
    val injection_number
    // tuple val(injection_number), val(tscrunch)
    // val tscrunch

    output:
    val injection_number

    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel

    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
    python3.6 ${projectDir}/pipeline_peasoup.py --tscrunch_index=2 --search_args=${params.search_params} --injection_file=${params.injection_plan} --out_dir=${params.output_dir}  --injection_number=${injection_number}
    inputs=""

    """
}

process peasoup3 {
    label "peasoup3"
    container "${params.peasoup_image}"

    input:
    val injection_number
    // tuple val(injection_number), val(tscrunch)
    // val tscrunch

    output:
    val injection_number

    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel

    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
    python3.6 ${projectDir}/pipeline_peasoup.py --tscrunch_index=3 --search_args=${params.search_params} --injection_file=${params.injection_plan} --out_dir=${params.output_dir}  --injection_number=${injection_number}
    inputs=""

    """
}

process candidate_filter {
    label "candidate_filter"
    container "${params.candidate_filter_image}"

    input:
    val injection_number

    output:
    val injection_number
    
    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel

    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
    python3.6 ${projectDir}/pipeline_candidate_filter.py --search_args=${params.search_params} --data_dir=${params.data_dir} --out_dir=${params.output_dir} --injection_number=${injection_number}
    
    """
}

process fold_cand {
    label "fold_cand"
    container "${params.fold_cand_image}"

    input:
    val injection_number

    output:
    val injection_number

    scratch "/tmp/rsenzel"

    script:
    """
    mkdir -p /tmp/rsenzel

    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/software/PulsarX/src/ymw16/.libs;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-10.2/targets/x86_64-linux/lib;
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib;
    python3.6 ${projectDir}/pipeline_fold.py --mode=cand --search_args=${params.search_params}  --injection_file=${params.injection_plan} --out_dir=${params.output_dir}  --injection_number=${injection_number} --ncpus=${task.cpus}

    """
}