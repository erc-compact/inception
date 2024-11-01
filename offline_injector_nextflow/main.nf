#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process injector{
    label "injector"
    container "${params.injector_image}"

    scratch true

    output:
    tuple path("candidates/*.candfile"), path("candidates/*_meta.txt")

    script:
    """
    #!/bin/bash
    # Ensure the output directory exists
    mkdir -p ${output_path}/Band_${params.band}/candidates

    # Running the XML Parsing task
    python3 ${params.injector_dir}/SCRIPT_inject_pulsar.py --signal ${params.inject_file} --fb {} --ncpu {} --ephem ${params.ephem_file} --output $get_cwd
    """
}



workflow {
    
}