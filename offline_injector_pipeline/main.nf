process pr1 {
    label "process_1"
    container "${params.sing_img}"
    publishDir "${params.pr1.output_path}", pattern: "*.csv", mode: 'copy'

    input:
    val number

    output:
    path "*.csv" 

    scratch "${params.pr1.tmp_dir}"

    script:
    """
    python3.6 ${projectDir}/test.py  --file ${params.pr1.inputs} --num $number
    """
}


workflow{
    numbers = Channel.from(0..99)
    pr1_out = pr1(numbers)

    pr1_out.view()
}