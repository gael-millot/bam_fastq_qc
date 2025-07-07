#!/usr/bin/env bash

#########################################################################
##                                                                     ##
##     Nremove.sh                                                      ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################



# $1 in_path
# $2 out_path
input_file=$1
output_file=$2
log=$3

echo -e "\nRemoval for reads made of N only\n" >> ${log}
awk '{lineKind=(NR-1)%4;}lineKind==0{record=$0; next}lineKind==1{toGet=!($0~/^N*$/); if(toGet) print record}toGet' ${input_file} > ${output_file}
# warning: with no output dir for log.txt, the file is created in \\wsl$\Ubuntu-20.04\home\gael\work\35\b826898b7be994ff13b7bc73bc88d8\
# get the bad sequences + 3 other lines of the fastq #see https://stackoverflow.com/questions/11793942/delete-lines-before-and-after-a-match-in-bash-with-sed-or-awk
# BEWARE: !/^(N*)$/ does not work to take the good seq, because the + line will be a good one and will print the 4 corresponding lines

LC_NUMERIC="en_US.UTF-8" # this is to have printf working for comma thousand separator
line_nb_before=$(cat ${input_file} | wc -l)
line_nb_after=$(cat ${output_file} | wc -l)
echo -e "\n\nNumber of sequences before removing reads with only N: $(printf "%'d" $((${line_nb_before} / 4)))\n" >> ${log}
echo -e "Number of sequences after removing reads with only N: $(printf "%'d" $((${line_nb_after} / 4)))\n" >> ${log}
echo -e "Ratio: " >> ${log}
echo -e $(printf '%.2f\n' $(echo $" $line_nb_after / $line_nb_before " | bc -l)) >> ${log} # the number in '%.2f\n' is the number of decimals
echo -e "\n\n" >> ${log}




