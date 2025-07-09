nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     main.nf                                                         ##
##     Comparative analysis of methylation in E. coli                  ##
##         using Hifi PacBio Revio long reads                          ##
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


process bam2fastq{
    label "bioconvert"

    input:
    path bam_ch // parallelization expected

    output:
    path "*.fastq", emit: fastq_ch
    path "*.log", emit: bam2fastq_log_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    FILENAME_INI=\$(basename -- "${bam_ch}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    echo -e "\n\n================\n\n\${FILENAME}\n\n================\n\n" > bam2fastq.log
    bioconvert \${FILENAME}.bam \${FILENAME}.fastq |& tee -a bam2fastq.log
    """
}


process Nremove { // remove the reads made of N only. See section 8.3 of the labbook 20200520
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
    cache 'true'

    input:
    path fastq_Nremove_ch // parallelization expected
    path kraken_db

    output:
    path(['*.kraken2', 'NULL']), emit: kraken_ch, optional: true
    path "kraken.log", emit: kraken_log_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    echo -e "\n\n================\n\n${fastq_Nremove_ch.baseName}\n\n================\n\n" > kraken.log
    if [[ ${system_exec} != "local" ]] ; then
        kraken2 --db ${kraken_db} --threads ${task.cpus} --report kraken.log > ${fastq_Nremove_ch.baseName}.kraken2
    else
        echo -e "\nNo kraken analysis performed in local running\n" > kraken.log
        echo "" > NULL
    fi
    """
}


process fastqc1 { // section 8.5 of the labbook 20200520
    label 'fastqc'
    publishDir "${out_path}/fastQC1", mode: 'copy', pattern: "*_fastqc.*", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path fastq_Nremove_ch

    output:
    path "fastqc1.log", emit: fastqc1_log_ch, optional: true
    path "*_fastqc.*", emit: fastqc1_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    echo -e "\n\n================\n\n${fastq_Nremove_ch.baseName}\n\n================\n\n" > fastqc1.log
    fastqc ${fastq_Nremove_ch} |& tee -a fastqc1.log
    """
}


process multiQC{
    label "multiqc"
    publishDir "${out_path}/reports", mode: 'copy', pattern: "multiqc_report.html", overwrite: false
    publishDir "${out_path}/reports", mode: 'copy', pattern: "multiqc.log", overwrite: false

    input:
    path fastq_Nremove_ch // no parallelization expected

    output:
    path "multiqc_report.html", emit: multiqc_ch

    script:
    """
    multiqc . -n multiqc_report.html |& tee -a multiqc.log
    """
}

process Q20 { // section 24.2 of the labbook 20200707
    label 'samtools' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', pattern: "q20.txt", overwrite: false
    // publishDir "${out_path}/files", mode: 'copy', pattern: "${file_name}_q20.bam", overwrite: false // 
    cache 'true'

    input:
    val file_name
    path bam from bowtie2_ch1

    output:
    path "${file_name}_q20_dup.bam", emit: q20_ch1, q20_ch2, q20_ch3, q20_ch4
    path "read_nb_before", emit: bow_read_nb_ch
    path "read_nb_after", emit: q20_read_nb_ch
    path "q20.txt"
    path "report.rmd", emit: log_ch11

    script:
    """
    samtools view -q 20 -b ${bam} > ${file_name}_q20_dup.bam |& tee q20.txt
    samtools index ${file_name}_q20_dup.bam
    echo -e "\\n\\n<br /><br />\\n\\n###  Q20 filtering\\n\\n" > report.rmd
    read_nb_before=\$(samtools view ${bam} | wc -l | cut -f1 -d' ') # -h to add the header
    read_nb_after=\$(samtools view ${file_name}_q20_dup.bam | wc -l | cut -f1 -d' ') # -h to add the header
    echo -e "\\n\\nNumber of sequences before Q20 filtering: \$(printf "%'d" \${read_nb_before})\\n" >> report.rmd
    echo -e "\\n\\nNumber of sequences after Q20 filtering: \$(printf "%'d" \${read_nb_after})\\n" >> report.rmd
    echo -e "Ratio: " >> report.rmd
    echo -e \$(printf "%.2f\n" \$(echo \$" \$read_nb_after / \$read_nb_before " | bc -l)) >> report.rmd # the number in '%.2f\\n' is the number of decimals
    echo -e "\\n\\n" >> report.rmd
    echo \$read_nb_before > read_nb_before # because nf cannot output values easily
    echo \$read_nb_after > read_nb_after
    """
}

process no_soft_clipping { // section 24.4 of the labbook 20200707
    label 'samtools' // see the withLabel: bash in the nextflow config file 
    cache 'true'

    input:
    path bam from q20_ch1

    output:
    path "report.rmd", emit: log_ch12

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n###  Control that no more soft clipping in reads\\n\\n" > report.rmd
    echo -e "nb of reads with soft clipping (S) in CIGAR: \$(printf "%'d" \$(samtools view ${bam} | awk '\$6 ~ /.*[S].*/{print \$0}' | wc -l | cut -f1 -d' '))" >> report.rmd
    echo -e "\\n\\ntotal nb of reads: \$(printf "%'d" \$(samtools view ${bam} | wc -l | cut -f1 -d' '))" >> report.rmd
    """
}

process duplicate_removal { // section 24.5 of the labbook 20200707. Warning: USING 5' AND 3' COORDINATES
    label 'samtools' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', pattern: "dup_report.txt", overwrite: false
    // publishDir "${out_path}/files", mode: 'copy', pattern: "${file_name}_q20_nodup.bam", overwrite: false // 
    cache 'true'

    input:
    val file_name
    path bam from q20_ch2
    path ref

    output:
    path "${file_name}_q20_nodup.bam", emit: dup_ch1, dup_ch2
    path "dup_read_nb", emit: dup_read_nb_ch
    path "dup_report.txt"
    path "report.rmd", emit: log_ch13

    script:
    """
    duplicate_removal.sh ${bam} ${ref} "${file_name}_q20_nodup.bam" "dup_report.txt" "report.rmd"
    """
}

process plot_read_length { // section 8.8 section 8.12 of the labbook 20200520
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/figures", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{plot_read_length_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path length2 from length_fastq_ini_ch
    path length3 from length_fastq_5p_filter_ch
    path length4 from length_fastq_5p_filter_cut_ch
    path length5 from length_cutoff_ch
    path cute_file

    output:
    path "*.png", emit: fig_ch2 // warning: several files
    path "plot_read_length_report.txt"
    path "report.rmd", emit: log_ch7

    script:
    """
    echo -e '
\\n\\n<br /><br />\\n\\n###  Length of initial reads\\n\\n
\\n\\n</center>\\n\\n
![Figure 2: Frequency of reads according to read size (in bp).](./figures/plot_read_length_ini.png){width=600}
\\n\\n</center>\\n\\n
\\n\\n<br /><br />\\n\\n###  Length of reads after selection of attC in 5 prime \\n\\n
\\n\\n</center>\\n\\n
![Figure 3: Frequency of reads according to read size (in bp).](./figures/plot_read_length_fivep_filtering.png){width=600}
\\n\\n</center>\\n\\n
\\n\\n<br /><br />\\n\\n###  Length of reads after trimming \\n\\n
\\n\\n</center>\\n\\n
![Figure 4: Frequency of reads according to read size (in bp).](./figures/plot_read_length_fivep_filtering_cut.png){width=600}
\\n\\n</center>\\n\\n
\\n\\n<br /><br />\\n\\n###  Read length after cut-off\\n\\n
\\n\\n</center>\\n\\n
![Figure 5: Frequency of reads according to read size (in bp).](./figures/plot_read_length_cutoff.png){width=600}
\\n\\n</center>\\n\\n
    ' > report.rmd
    plot_read_length.R "${length2}" "${length3}" "${length4}" "${length5}" "${cute_file}" "plot_read_length_report.txt"
    """
    // single quotes required because of the !
}




process coverage { // section 24.5 of the labbook 20200707. Warning: USING 5' AND 3' COORDINATES
    label 'bedtools' // see the withLabel: bash in the nextflow config file 
    // publishDir "${out_path}/reports", mode: 'copy', pattern: "cov_report.txt", overwrite: false // inactivated because no cov_report published in "${out_path}/reports" probably because of the parallelization
    publishDir "${out_path}/files", mode: 'copy', pattern: "*.cov", overwrite: false
    cache 'true'

    input:
    path bam from bowtie2_ch2.concat(q20_ch3, dup_ch1)
    // file ref from ref_ch3 // not required because bedtools genomecov-g ${ref} not required when inputs are bam files

    output:
    path "*_mini.cov", emit: cov_ch // warning: 3 files
    // file "*.cov" // coverage per base if ever required but long process
    path "cov_report.txt", emit: cov_report_ch

    script:
    """
    # bedtools genomecov -d -ibam \${bam} > \${bam.baseName}.cov |& tee cov_report.txt # coverage per base if ever required but long process
    # to add the chr names | awk '{h[\$NF]++}; END { for(k in h) print k, h[k] }' | sort -V > \${bam.baseName}.cov
    bedtools genomecov -bga -ibam ${bam}  > ${bam.baseName}_mini.cov |& tee cov_report.txt
    # -g \${ref} not required when inputs are bam files
    """
}


//cov_report_ch.collectFile(name: "cov_report.txt").subscribe{it -> it.copyTo("${out_path}/reports")} // concatenate all the cov_report.txt files in channel cov_report_ch into a single file published into 



process plot_coverage { // section 24.6 of the labbook 20200707
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/figures", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    // publishDir "${out_path}/reports", mode: 'copy', pattern: "{plot_coverage_report.txt}", overwrite: false // 
    cache 'true'

    input:
    val file_name
    path cov from cov_ch // warning: 3 files
    path read_nb from bow_read_nb_ch.concat(q20_read_nb_ch, dup_read_nb_ch)
    val ori_coord
    val ter_coord
    val color_coverage
    val xlab
    path cute_file

    output:
    path "plot_${cov.baseName}.png", emit: fig_ch3 // warning: several files
    path "plot_coverage_report.txt", emit: plot_cov_report_ch

    script:
    """
    plot_coverage.R "${cov.baseName}" "${read_nb}" "${ori_coord}" "${ter_coord}" "${color_coverage}" "${xlab}" "${file_name}" "${cute_file}" "plot_coverage_report.txt"
    """
    // single quotes required because of the !
}


// plot_cov_report_ch.collectFile(name: "plot_cov_report.txt").subscribe{it -> it.copyTo("${out_path}/reports")} // concatenate all the cov_report.txt files in channel cov_report_ch into a single file published into ${out_path}/reports


process print_report{
    label 'r_ext'

    publishDir path: "${out_path}", mode: 'copy', pattern: "{*.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{print_report.log}", overwrite: false
    cache 'false'

    input:
    path template_rmd
    path ini_read_nb
    path n_remove_read_nb
    val nb_bam_files

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
            nb_bam_files = ${nb_bam_files}
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
    echo -e "full .nextflow.log is in: ${launchDir}\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
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
    if( ! (file("${projectDir}/bin/coverage.sh").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE coverage.sh FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    // AB_model not trested because in parameters
    if( ! (file("${projectDir}/bin/plot_coverage.R").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE plot_coverage.R FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (file("${projectDir}/bin/html_report_template.rmd").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE html_report_template.rmd FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
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
            tempo_log = files.every { it.name.toLowerCase().endsWith('.bam') }
            if ( ! tempo_log) {
                error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE.\nTHE INDICATED FOLDER CAN ONLY CONTAIN .bam FILES: ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
            }
        }
    }else if( ! (sample_path =~ /.*\.(zip|bam)$/)){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE.\nIF NOT A FOLDER, MUST BE A FILE FINISHING BY .zip OR .bam.\nHERE IT IS: ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (ref_path in String || ref_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ref_path PARAMETER IN nextflow.config FILE:\n${ref_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(ref_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ref_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${ref_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if("${system_exec}" != "local"){
        if( ! (kraken_db_path in String || kraken_db_path in GString) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID kraken_db_path PARAMETER IN nextflow.config FILE:\n${kraken_db_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
        }else if( ! (file(kraken_db_path).exists()) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID kraken_db_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${kraken_db_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }
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

    file("${out_path}/fastQC1").mkdirs()
    file("${out_path}/files").mkdirs()

    //////// end Folder creation

    //////// files import

    cute_file = file(cute_path) // in variable because a single file
    template_rmd = file(template_rmd_path)
    if(system_exec != 'local'){
        kraken_db = file(kraken_db_path)
    }else{
        kraken_db = file("NULL_kraken_db")
    }


    //////// end files import


    //////// Main

    if(sample_path =~ /.*\.zip$/){
        Unzip( // warning: it is a process defined above
            Channel.fromPath(sample_path)
        ) 
        bam_ch = Unzip.out.unzip_ch.flatten()
    }else if(sample_path =~ /^http.*\.bam$/){
        bam_ch = Channel.fromPath("${sample_path}", checkIfExists: false) // in channel because many files 
    }else{
        bam_ch = Channel.fromPath("${sample_path}/*.*", checkIfExists: false) // in channel because many files 
    }
    
    nb_bam_files = bam_ch.count()

    workflowParam(
        modules
    )

    bam2fastq(
        bam_ch
    )
    bam2fastq.out.bam2fastq_log_ch.collectFile(name: "bam2fastq.log").subscribe{it -> it.copyTo("${out_path}/reports")}

    Nremove(
        bam2fastq.out.fastq_ch
    )
    ini_read_nb_ch2 = Nremove.out.ini_read_nb_ch.collectFile(name: "ini_read_nb.tsv", skip: 1, keepHeader: true)
    ini_read_nb_ch2.subscribe{it -> it.copyTo("${out_path}/files")}
    n_remove_read_nb_ch2 = Nremove.out.n_remove_read_nb_ch.collectFile(name: "n_remove_read_nb.tsv", skip: 1, keepHeader: true)
    n_remove_read_nb_ch2.subscribe{it -> it.copyTo("${out_path}/files")}

    kraken(
        Nremove.out.fastq_Nremove_ch,
        kraken_db
    )
    kraken.out.kraken_log_ch.collectFile(name: "kraken.log").subscribe{it -> it.copyTo("${out_path}/reports")}
    if(system_exec != 'local'){
        kraken.out.kraken_ch.collectFile(name: "kraken_report.html").subscribe{it -> it.copyTo("${out_path}/files")}
    }

    fastqc1(
        Nremove.out.fastq_Nremove_ch
    )
    fastqc1.out.fastqc1_log_ch.collectFile(name: "fastqc1.log").subscribe{it -> it.copyTo("${out_path}/reports")}

    multiQC(
        fastqc1.out.fastqc1_ch.mix(kraken.out.kraken_ch).collect()
    )

    print_report(
        template_rmd, 
        ini_read_nb_ch2, 
        n_remove_read_nb_ch2, 
        nb_bam_files
    )

// bam_ch.collect().ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE distance_hist PROCESS\n\n========\n\n"}

    backup(
        config_file, 
        log_file
    )

}
