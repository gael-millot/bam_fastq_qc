nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     main.nf                                                         ##
##     bam_fastq_qc                                                    ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################
*/


//////// Processes


process workflowParam { // create a file with the workflow parameters in out_path
    label 'bash'
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false
    cache 'false'

    input:
    val modules

    output:
    path "Run_info.txt"

    script:
    """
    echo "Project (empty means no .git folder where the main.nf file is present): " \$(git -C ${projectDir} remote -v | head -n 1) > Run_info.txt # works only if the main script run is located in a directory that has a .git folder, i.e., that is connected to a distant repo
    echo "Git info (empty means no .git folder where the main.nf file is present): " \$(git -C ${projectDir} describe --abbrev=10 --dirty --always --tags) >> Run_info.txt # idem. Provide the small commit number of the script and nextflow.config used in the execution
    echo "Cmd line: ${workflow.commandLine}" >> Run_info.txt
    echo "execution mode": ${system_exec} >> Run_info.txt
    modules=$modules # this is just to deal with variable interpretation during the creation of the .command.sh file by nextflow. See also \$modules below
    if [[ ! -z \$modules ]] ; then
        echo "loaded modules (according to specification by the user thanks to the --modules argument of main.nf): ${modules}" >> Run_info.txt
    fi
    echo "Manifest's pipeline version: ${workflow.manifest.version}" >> Run_info.txt
    echo "result path: ${out_path}" >> Run_info.txt
    echo "nextflow version: ${nextflow.version}" >> Run_info.txt
    echo -e "\\n\\nIMPLICIT VARIABLES:\\n\\nlaunchDir (directory where the workflow is run): ${launchDir}\\nprojectDir (directory where the main.nf script is located): ${projectDir}\\nworkDir (directory where tasks temporary files are created): ${workDir}" >> Run_info.txt
    echo -e "\\n\\nUSER VARIABLES:\\n\\nout_path: ${out_path}\\nsample_path: ${sample_path}" >> Run_info.txt
    """
}
//${projectDir} nextflow variable
//${workflow.commandLine} nextflow variable
//${workflow.manifest.version} nextflow variable
//Note that variables like ${out_path} are interpreted in the script block



process Unzip {
    label 'unzip'
    cache 'true'

    input:
    path zip

    output:
    path "*.bam", emit: unzip_ch

    script:
    """
    #!/bin/bash -ue
    FILENAME_INI=\$(basename -- "${zip}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    if [[ ! "\${FILE_EXTENSION}" =~ zip ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE FILE EXTENSION MUST BE \\".zip\\" AND NOT \${FILENAME_INI}\\n\\n========\\n\\n"
        exit 1
    fi
    TEMPO=\$(unzip -l ${zip} | awk 'NR>3 {print \$4}' | grep -v '^\$') # get the file list in the archive
    # Check if only files in the archives and no directory
    for file in "\$TEMPO" ; do
        if [[ ! "\$file" =~ \\.(bam)\$ ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nALL THE ELEMENTS IN THE .zip FILE MUST BE .bam FILES ONLY.\\nHERE AN ELEMENT IS:\\n\$file.\\n\\n========\\n\\n"
            exit 1
        fi
    done
    # end Check if only files in the archives and no directory
    unzip ${zip} 
    """
}


process bam2fastq{ // is indicated as bam2fastq - during nextflow running if bam_ch is empty
    label "bioconvert"

    input:
    path bam_ch // parallelization expected

    output:
    path "*.fastq", emit: tempo_fastq_ch
    path "*.log", emit: bam2fastq_log_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    FILENAME_INI=\$(basename -- "${bam_ch}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    echo -e "\\n\\n================\\n\\n\${FILENAME}\\n\\n================\\n\\n" > bam2fastq.log
    bioconvert \${FILENAME}.bam \${FILENAME}.fastq |& tee -a bam2fastq.log
    """
}

process fastqc {
    label 'fastqc'
    //publishDir "${out_path}/fastQC1", mode: 'copy', pattern: "*_fastqc.*", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path fastq // parallelization expected

    output:
    path "fastqc.log", emit: fastqc_log_ch, optional: true
    path "*_fastqc.*", emit: fastqc_print_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    echo -e "\\n\\n================\\n\\n${fastq.baseName}\\n\\n================\\n\\n" > fastqc.log
    fastqc ${fastq} |& tee -a fastqc.log
    """
}

process Nremove { // remove the reads made of N only.
    label 'bash' // see the withLabel: bash in the nextflow config file 
    cache 'true'

    input:
    path fastq_ch // parallelization expected

    output:
    path "*_Nremove.fastq", emit: fastq_Nremove_ch
    path "ini_read_nb.tsv", emit: ini_read_nb_ch
    path "n_remove_read_nb.tsv", emit: n_remove_read_nb_ch

    script:
    """
    #!/bin/bash -ue
    FILENAME_INI=\$(basename -- "${fastq_ch}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    awk '{lineKind=(NR-1)%4;}lineKind==0{record=\$0; next}lineKind==1{toGet=!(\$0~/^N*\$/); if(toGet) print record}toGet' ${fastq_ch} > \${FILENAME}_Nremove.fastq
    # get the bad sequences + 3 other lines of the fastq #see https://stackoverflow.com/questions/11793942/delete-lines-before-and-after-a-match-in-bash-with-sed-or-awk
    # BEWARE: !/^(N*)\$/ does not work to take the good seq, because the + line will be a good one and will print the 4 corresponding lines
    LC_NUMERIC="en_US.UTF-8" # this is to have printf working for comma thousand separator
    line_nb_before=\$(cat ${fastq_ch} | wc -l)
    line_nb_after=\$(cat \${FILENAME}_Nremove.fastq | wc -l)
    echo -e "FILE_NAME\\tREADS_NB" > ini_read_nb.tsv
    echo -e "\${FILENAME_INI}\\t\$(printf "%'d" \$((\${line_nb_before} / 4)))" >> ini_read_nb.tsv
    echo -e "FILE_NAME\\tREADS_NB\\tRATIO" > n_remove_read_nb.tsv
    echo -e "\${FILENAME_INI}\\t\$(printf "%'d" \$((\${line_nb_after} / 4)))\\t\$(printf '%.2f\\n' \$(echo \$" \$line_nb_after / \$line_nb_before " | bc -l))" >> n_remove_read_nb.tsv
    # the number in '%.2f\\n' is the number of decimals
    """
}


process kraken {
    label 'kraken'
    publishDir "${out_path}/kraken", mode: 'copy', pattern: "*_kraken2.txt", overwrite: false
    cache 'true'

    input:
    path fastq_Nremove_ch // parallelization expected
    path kraken_db

    output:
    path "*_kraken2.txt", emit: kraken_ch, optional: true
    path "kraken.log", emit: kraken_log_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    echo -e "\\n\\n================\\n\\n${fastq_Nremove_ch.baseName}\\n\\n================\\n\\n" > kraken.log
    free -h # display the memory available
    kraken2 --db ${kraken_db} --threads \$(nproc) --report ${fastq_Nremove_ch.baseName}_report_kraken2.txt ${fastq_Nremove_ch} > ${fastq_Nremove_ch.baseName}_classif_kraken2.txt 2> kraken.log # 2> kraken.log because if the warnings or errors are in \${fastq_Nremove_ch.baseName}_report_kraken2.txt, multiqc will not recognize the file
    """
}

process multiQC{
    label "multiqc"
    publishDir "${out_path}/multiQC", mode: 'copy', pattern: "multiqc_report*", overwrite: false
    publishDir "${out_path}/reports", mode: 'copy', pattern: "multiqc.log", overwrite: false

    input:
    path fastq_Nremove_ch // no parallelization expected

    output:
    path "multiqc_report*"
    path "multiqc.log"

    script:
    """
    multiqc . --filename multiqc_report.html |& tee -a multiqc.log
    """
}

process read_cutoff {
    label 'bash'
    cache 'true'

    input:
    path fastq_Nremove_ch // parallelization expected
    val cutoff

    output:
    path "*_cutoff.fastq", emit: fastq_cutoff_ch
    path "cutoff_read_nb.tsv", emit: cutoff_read_nb_ch


    script:
    """
    #!/bin/bash -ue
    FILENAME_INI=\$(basename -- "${fastq_Nremove_ch}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    # cutoff
    awk -v var1=${cutoff} '{lineKind=(NR-1)%4}lineKind==0{record=\$0; next}lineKind==1{toGet=(length(\$0)>=var1); if(toGet) print record}toGet' ${fastq_Nremove_ch} > \${FILENAME}_cutoff.fastq
    # end cutoff
    LC_NUMERIC="en_US.UTF-8" # this is to have printf working for comma thousand separator
    line_nb_before=\$(cat ${fastq_Nremove_ch} | wc -l)
    line_nb_after=\$(cat \${FILENAME}_cutoff.fastq | wc -l)
    echo -e "FILE_NAME\\tREADS_NB\\tRATIO" > cutoff_read_nb.tsv
    echo -e "\${FILENAME_INI}\\t\$(printf "%'d" \$((\${line_nb_after} / 4)))\\t\$(printf '%.2f\\n' \$(echo \$" \$line_nb_after / \$line_nb_before " | bc -l))" >> cutoff_read_nb.tsv
    # the number in '%.2f\\n' is the number of decimals
    """
}


process read_length {
    label 'bash'
    cache 'true'

    input:
    path fastq // parallelization expected

    output:
    path "*_read_length.tsv", emit: read_length_ch

    script:
    """
    #!/bin/bash -ue
    FILENAME_INI=\$(basename -- "${fastq}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    # length
    awk -v var1=\${FILENAME_INI} 'BEGIN{print "FILE_NAME\\tREAD_NAME\\tLENGTH"}{lineKind=(NR-1)%4}lineKind==0{record=\$0; next}lineKind==1{print var1"\\t"record"\\t"length(\$0)}' ${fastq} > \${FILENAME}_read_length.tsv
    # end length
    """
}


process plot_read_length {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/figures", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{plot_read_length_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/files", mode: 'copy', pattern: "{*_freq.tsv}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path tsv // parallelization expected
    path cute_file

    output:
    path "*.png", emit: png_ch
    path "*_freq.tsv"
    path "plot_read_length_report.txt"

    script:
    """
    plot_read_length.R  "${tsv}" "${cute_file}" "plot_read_length_report.txt"
    """
    // single quotes required because of the !
}

process print_report{
    label 'r_ext'

    publishDir path: "${out_path}", mode: 'copy', pattern: "{*.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{print_report.log}", overwrite: false
    cache 'false'

    input:
    path template_rmd
    path ini_read_nb
    path n_remove_read_nb
    path cutoff_read_nb_ch2
    path png_ch
    val nb_bam_files
    val cutoff

    output:
    path "report.html"
    path "print_report.log"

    script:
    """
    #!/bin/bash -ue
    cp ${template_rmd} report_file.rmd
    cp -r "${out_path}/files" .
    Rscript -e '
        rmarkdown::render(
        input = "report_file.rmd",
        output_file = "report.html",
        # list of the variables waiting to be replaced in the rmd file :
        params = list(
            nb_bam_files = ${nb_bam_files}, 
            cutoff = ${cutoff}
        ),
        # output_dir = ".",
        # intermediates_dir = "./",
        # knit_root_dir = "./",
        run_pandoc = TRUE,
        quiet = TRUE,
        clean = TRUE
        )
    ' |& tee -a print_report.log
    """
}

// Save the config file and the log file for a specific run
process backup {
    label 'bash'
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false // since I am in mode copy, all the output files will be copied into the publishDir. See \\wsl$\Ubuntu-20.04\home\gael\work\aa\a0e9a739acae026fb205bc3fc21f9b
    cache 'false'

    input:
    path config_file
    path log_file

    output:
    path "${config_file}" // warning message if we use path config_file
    path "${log_file}" // warning message if we use path log_file
    path "Log_info.txt"

    script:
    """
    #!/bin/bash -ue
    echo -e "full .nextflow.log is in: ${launchDir}\\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}

//////// End Processes

//////// Workflow



workflow {

    print("\n\nINITIATION TIME: ${workflow.start}")

    //////// Options of nextflow run

    // --modules (it is just for the process workflowParam)
    params.modules = "" // if --module is used, this default value will be overridden
    // end --modules (it is just for the process workflowParam)

    //////// end Options of nextflow run


    //////// Variables

    modules = params.modules // remove the dot -> can be used in bash scripts
    config_file = workflow.configFiles[0] // better to use this than config_file = file("${projectDir}/nextflow.config") because the latter is not good if -c option of nextflow run is used
    log_file = file("${launchDir}/.nextflow.log")

    //////// end Variables


    //////// Checks
    //// check of the bin folder
    if( ! (file("${projectDir}/bin/html_report_template.rmd").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE html_report_template.rmd FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    // AB_model not trested because in parameters
    if( ! (file("${projectDir}/bin/plot_read_length.R").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE plot_read_length.R FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (file("${projectDir}/bin/cute_little_R_functions_12.8.R").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE cute_little_R_functions_12.8.R FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    //// end check of the bin folder
    if( ! (sample_path in String || sample_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE:\n${sample_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(sample_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }else if (file(sample_path).isDirectory()) {
        files = file(sample_path).listFiles().findAll { it.isFile() }
        if (files.isEmpty()) {
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE.\nCANNOT BE AN EMPTY FOLDER: ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }else{
            tempo_log = files.every{it.name.toLowerCase() ==~ /.*\.(bam|fq|fastq)$/} // for single end: files.every{it.name.toLowerCase().endsWith('.bam') }
            if ( ! tempo_log) {
                error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE.\nTHE INDICATED FOLDER CAN ONLY CONTAIN .bam FILES: ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
            }
        }
    }else if( ! (sample_path =~ /.*\.(zip|bam|fq|fastq)$/)){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE.\nIF NOT A FOLDER, MUST BE A FILE FINISHING BY .zip OR .bam OR .fq OR .fastq.\nHERE IT IS: ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (kraken_db_path in String || kraken_db_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID kraken_db_path PARAMETER IN nextflow.config FILE:\n${kraken_db_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(kraken_db_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID kraken_db_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${kraken_db_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (cutoff in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cutoff PARAMETER IN nextflow.config FILE:\n${cutoff}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (cutoff =~  /^[0-9]+$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cutoff PARAMETER IN nextflow.config FILE:\n${cutoff}\nMUST BE A POSITIVE INGETER\n\n========\n\n"
    }
    if( ! (cute_path in String || cute_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN nextflow.config FILE:\n${cute_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(cute_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${cute_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }

    // below : those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
    // out_ini
    print("\n\nRESULT DIRECTORY: ${out_path}")
    //print("\n\nWARNING: PARAMETERS ALREADY INTERPRETED IN THE .config FILE:")
    //print("    system_exec: ${system_exec}")
    //print("    out_path: ${out_path_ini}")
    if("${system_exec}" == "slurm"){
        print("    queue: ${slurm_queue}")
        print("    qos: ${slurm_qos}")
    }
    if("${system_exec}" != "local"){
        print("    add_options: ${add_options}")
    }
    print("\n\n")


    //////// end Checks


    //////// Variable modification


    //////// end Variable modification


    //////// Channels

    // bam_ch define below because can be a .zip file

    //////// end Channels

    //////// Folder creation

    file("${out_path}/fastQC").mkdirs()
    file("${out_path}/multiQC").mkdirs()
    file("${out_path}/kraken").mkdirs()
    file("${out_path}/files").mkdirs()
    file("${out_path}/figures").mkdirs()

    //////// end Folder creation

    //////// files import

    cute_file = file(cute_path) // in variable because a single file
    template_rmd = file(template_rmd_path)
    kraken_db = file(kraken_db_path)


    //////// end files import


    //////// Main

    workflowParam(
        modules
    )

    if(sample_path =~ /.*\.zip$/){
        Unzip( // warning: it is a process defined above
            Channel.fromPath(sample_path)
        ) 
        files_ch = Unzip.out.unzip_ch.flatten()
    }else if(sample_path =~ /^http.*\.(bam|fq|fastq)$/){
        files_ch = Channel.fromPath("${sample_path}", checkIfExists: false) // in channel because many files 
    }else{
        files_ch = Channel.fromPath("${sample_path}/*.*", checkIfExists: false) // in channel because many files 
    }
    
    nb_bam_files = files_ch.count()

    bam_ch = files_ch.filter { it.name ==~ /.*\.bam(\.gz)?$/ }
    ini_fastq_ch = files_ch.filter { it.name ==~ /.*\.(fastq|fq)(\.gz)?$/ }

    bam2fastq( // is indicated as bam2fastq - during nextflow running if bam_ch is empty
        bam_ch
    )
    bam2fastq.out.bam2fastq_log_ch.collectFile(name: "bam2fastq.log").subscribe{it -> it.copyTo("${out_path}/reports")}

    fastq_ch = ini_fastq_ch.mix(bam2fastq.out.tempo_fastq_ch)

    fastqc(
        fastq_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE bam2fastq PROCESS\n\n========\n\n"}
    )
    fastqc.out.fastqc_log_ch.collectFile(name: "fastqc1.log").subscribe{it -> it.copyTo("${out_path}/reports")}
    fastqc.out.fastqc_print_ch.flatten().subscribe{it -> it.copyTo("${out_path}/fastQC")}


    Nremove(
        fastq_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE bam2fastq PROCESS\n\n========\n\n"}
    )
    ini_read_nb_ch2 = Nremove.out.ini_read_nb_ch.collectFile(name: "ini_read_nb.tsv", skip: 1, keepHeader: true)
    ini_read_nb_ch2.subscribe{it -> it.copyTo("${out_path}/files")}
    n_remove_read_nb_ch2 = Nremove.out.n_remove_read_nb_ch.collectFile(name: "n_remove_read_nb.tsv", skip: 1, keepHeader: true)
    n_remove_read_nb_ch2.subscribe{it -> it.copyTo("${out_path}/files")}


    kraken(
        Nremove.out.fastq_Nremove_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Nremove PROCESS\n\n========\n\n"},
        kraken_db
    )
    kraken.out.kraken_log_ch.collectFile(name: "kraken.log").subscribe{it -> it.copyTo("${out_path}/reports")}
 
    multiQC(
        fastqc.out.fastqc_print_ch.mix(kraken.out.kraken_ch).collect().ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE kraken AND fastqc PROCESSES\n\n========\n\n"}
    )


    read_cutoff(
        Nremove.out.fastq_Nremove_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Nremove PROCESS\n\n========\n\n"},
        cutoff
    )
    cutoff_read_nb_ch2 = read_cutoff.out.cutoff_read_nb_ch.collectFile(name: "cutoff_read_nb.tsv", skip: 1, keepHeader: true)
    cutoff_read_nb_ch2.subscribe{it -> it.copyTo("${out_path}/files")}


    read_length(
        Nremove.out.fastq_Nremove_ch.mix(read_cutoff.out.fastq_cutoff_ch).ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Nremove AND read_cutoff PROCESSES\n\n========\n\n"}
    )
    read_length_ch2 = read_length.out.read_length_ch.collectFile(name: "read_length.tsv", skip: 1, keepHeader: true)
    read_length_ch2.subscribe{it -> it.copyTo("${out_path}/files")}


    plot_read_length(
        read_length_ch2,
        cute_file
    )


    print_report(
        template_rmd, 
        ini_read_nb_ch2, 
        n_remove_read_nb_ch2, 
        cutoff_read_nb_ch2, 
        plot_read_length.out.png_ch, 
        nb_bam_files,
        cutoff
    )


    backup(
        config_file, 
        log_file
    )

}
