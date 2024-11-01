#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process parse_xml{
    label "xml_parse"
    container "${params.pulsarx_image}"
    publishDir "${params.output_path}/Band_${params.band}/", pattern: "**/*.candfile", mode: 'copy'
    publishDir "${params.output_path}/Band_${params.band}/", pattern: "**/*_meta.txt", mode: 'copy'

    input:
    path xml_file
    val output_path


    output:
    tuple path("candidates/*.candfile"), path("candidates/*_meta.txt")

    script:
    """
    #!/bin/bash
    # Ensure the output directory exists
    mkdir -p ${output_path}/Band_${params.band}/candidates

    # Running the XML Parsing task
    python3 ${baseDir}/candidate_parser.py -i ${xml_file} -o candidates/  -f ${params.fold_technique} -cpn ${params.pulsarx_cpn} -b ${params.beam_name} -rfi ${params.rfi_filter} -n ${params.nh} -nsub ${params.nsub} -clfd ${params.clfd} -fn ${params.fast_nbins} -sn ${params.slow_nbins}
    """
}


process pulsarx_fold{
    label "pulsarx_fold"
    container "${params.pulsarx_image}"
    publishDir "${params.output_path}/Band_${params.band}/", pattern: "${params.beam_name}_*/**", mode: 'copy'
    
    input:
    path xml_parser_results
    val output_path

    output:
    tuple path("**/*.png"), path("**/*.ar"), path("**/*.cands")

    script:
    """
    #!/bin/bash
    python3 ${baseDir}/pulsarx_fold.py -threads ${params.threads} -p ${params.template_dir} -meta ${xml_parser_results[1]} -cands ${xml_parser_results[0]} -fits ${params.isitfits} -band ${params.band}
    """
}

process filtool {
    label "filtool"
    container "${params.pulsarx_image}"
    publishDir "${params.filtool.output_path}", pattern: "*.fil", mode: 'copy' 
    
    input:
    path file_path

    output:
    path("*.fil")

    when:
    def file_path_str = file_path.toString()
    // check the file name and see if _Band_ is present, don't proceed for that file
    file_path_str.contains("_Band_")
    script:
    def file_path_str = file_path.toString()
    def file_name = new File(file_path_str).name.replaceAll(/\.fits$|\.sf$|\.rf$/, "")
    def GC = new File(file_path_str).name.replaceAll(/_Band_.*$/, "")
    def Band = new File(file_path_str).name.replaceAll(/.*_Band_(\d+).*/, '$1')
    def mask_param_name = "Band${Band}_Mask"
    // println "mask_param_name: $mask_param_name
    // Get the value of the mask using the dynamic mask_param_name
    def mask_value = params.filtool[mask_param_name]
    if (mask_value == null) {
        throw new IllegalStateException("Mask parameter ${mask_param_name} is not defined in params.filtool")
    }
    // println "mask_value: $mask_value"

    """
    #!/bin/bash
    # Run filtool with extracted values
    filtool -v --flip -o ${file_name} --psrfits ${file_path} -z ${mask_value}
    """
}

process peasoup {
    label "peasoup"
    container "${params.peasoup_image}"
    publishDir "${params.peasoup.output_path}", pattern: "**/*.xml", mode: 'copy'
    
    publishDir "${params.peasoup.output_path}", pattern: "**/*.xml", mode: 'copy'
    

    input:
    tuple val(fil_file), val(fft_size)
    each path(dm_file)

    output:
    tuple path("**/*.xml") , env(dm_name)

    script:
    def file_path_str = fil_file.toString()
    def GC = new File(file_path_str).name.replaceAll(/_Band_.*$/, "")
    def bandName = new File(file_path_str).name.replaceAll(/.*_Band_(\d+).*/, '$1')
    def bandName = new File(file_path_str).name.replaceAll(/.*_Band_(\d+).*/, '$1')
    """
    #!/bin/bash
    dm_name=\$(basename "${dm_file}" .dm)

    echo "Running peasoup with DM file: \${dm_name}" and FFT size: ${fft_size} on ${fil_file}

    peasoup -p -i ${fil_file} -o ${GC}_${bandName}_\${dm_name} --limit 100000 -m ${params.peasoup.min_snr} -t 1 --acc_start ${params.peasoup.acc_start} --acc_end ${params.peasoup.acc_end} --dm_file ${dm_file} --ram_limit_gb ${params.peasoup.ram_limit_gb} -n ${params.peasoup.nharmonics} --fft_size ${fft_size}
    """
    
}



process temp{
    script1:
    """
    filtool -v -t {num_threads} --zapthre {zapping_threshold} --fd {fscrunch} --filplan {filplan_fname} -l {segment_length} 
            --baseline {0} {0} --fillPatch rand -z {rfi_flags} -o {rootname} -f {' '.join(input_fils)}"
    """

    script2:
    """
    peasoup -k {channel_mask} -z {birdie_list} "
           f"-i {input_fil} --ram_limit_gb {ram_limit} "
           f"--dm_file {dm_list} --limit {candidate_limit} "
           f"-n {nharmonics}  -m {snr_threshold} --acc_start {start_accel} "
           f"--acc_end {end_accel} --fft_size {fft_length} -o {out_dir} --dedisp_gulp {gulp_size}
    """

    script3:
    """
    candidate_filter.py -i %s -o %s/%d --threshold %f -c /home/psr/software/candidate_filter/candidate_filter/default_config.json 
                        --rfi /home/psr/software/candidate_filter/candidate_filter/known_rfi.txt --p_tol %f --dm_tol %f" %
                        (candidate_files_list_path, tmp_dir, processing_id, snr_cutoff, processing_args['p_tol'], processing_args['dm_tol']), shell=True)
    """

    script4:
    """
    psrfold_fil2 --dmboost 250 --plotx -v -t 12 --candfile {} -n {} {} {} --template {} --clfd 8 -L {} --fillPatch rand -f {} --rfi zdot {} --fd {} --td {}".format(
                    pred_file, nsubband, nbins_string, beam_tag, TEMPLATE, subint_length, input_filenames, zap_string, fscrunch, tscrunch)
    """
}

workflow {
    xml_queue = Channel.fromPath(params.xml_files)
    xml_parser_results= parse_xml(xml_queue, params.output_path).transpose()
    
    pulsarx_fold(xml_parser_results, params.output_path)
}

workflow {
    fits_files_ch = find_files()
    fits_queue = fits_files_ch.splitText().map { it.trim() }
    fits_queue.view()
    filtool_channel = filtool(fits_queue)
    // filtool_channel = Channel.fromPath("/hercules/scratch/fkareem/NGC7099/Filtool/NGC7099_Band_*fil") //This is a quick fix to avoid the filtool process for fil files
    nearest_power_two_results = nearest_power_of_two_calculator(filtool_channel).transpose()
    filtool_output = nearest_power_two_results.map { item -> 
        def (fil_file, fft_size) = item
        return [fil_file, fft_size]
        }
    dmfiles = dmfile_gen()
    peasoup_results = peasoup(filtool_output, dmfiles)
    xml_parser_results= parse_xml(peasoup_results).transpose()
    xml_parser_results.view()
    fold_out = pulsarx_fold(xml_parser_results).collectFile(name: 'finished.txt', newLine : true)
    pics_classifier(fold_out)
}
