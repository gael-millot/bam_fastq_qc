#!/usr/bin/env bash

#########################################################################
##                                                                     ##
##     coverage.sh                                                     ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################



# $1 in_path
# $2 out_path
bam1=$1
bam2=$2
bam3=$3
ref=$4
file_name=$5
exec_log=$6
log=$7

# bam=test.fastq2_bowtie2_sorted_q20.bam
# ref_genome=Ecoli-K12-MG1655_ORI_CENTERED.fasta
# output_file="test.fastq2_nodup_sorted.bam"
# exec_log="cov_report.txt"
# log="report.rmd"


echo -e "<br /><br />\n\n### Coverage computation\n\n" >> ${log}


bedtools genomecov -d -ibam ${bam1}  -g ${ref} > ${file_name}.cov
# to add the chr names | awk '{h[$NF]++}; END { for(k in h) print k, h[k] }' | sort -V > ${file_name}.cov.txt
bedtools genomecov -d -ibam ${bam2}  -g ${ref} > ${file_name}_q20.cov

bedtools genomecov -d -ibam ${bam3}  -g ${ref} > ${file_name}_q20_nodup.cov


bedtools genomecov -bga -ibam ${file_name}.bam  -g ${ref} > ${file_name}_mini.cov
# to add the chr names | awk '{h[$NF]++}; END { for(k in h) print k, h[k] }' | sort -V > ${file_name}.mini.cov.txt
bedtools genomecov -bga -ibam ${file_name}_q20.bam  -g ${ref} > ${file_name}_q20_mini.cov

bedtools genomecov -bga -ibam ${file_name}_q20.nodup.sorted.bam  -g ${ref} > ${file_name}_q20_nodup_mini.cov


















echo -e "\n\n################ Check that no BX:Z: TAG already exists\n\n" >> ${exec_log}
samtools view ${bam} | awk 'BEGIN{TEST=1}/^@/{next}/^.*BX:Z:.*$/{print "\n\n============\n\nERROR: BX:Z: TAG ALREADY EXISTS\n\n============\n\n" ; TEST=0 ; exit 1}END{if(TEST==1){print "\n\nNO BX:Z: TAG PRESENT\n\n"}}' |& tee -a ${exec_log}

# creates THE BX:Z: TAG which is FLAG_POS_(POS + length(SEQ) - 1)
echo -e "\n\n################ Creates THE BX:Z: TAG which is FLAG_POS_(POS + length(SEQ) - 1)\n\n" >> ${exec_log}
samtools view -h ${bam} | awk -v var1=4 -v var2=$(awk '{lineKind=(NR-1)%2}lineKind==1{print length(length($0))}' $ref_genome) 'BEGIN{FS="\t" ; OFS="" ; ORS=""}/^@/{print $0"\n" ; next}{print $0"\tBX:Z:" ; printf("%0"var1"d", $2) ; print "_" ; printf("%0"var2"d", $4) ; print "_" ; printf("%0"var2"d", $4 + length($10) - 1) ; print "\n"}' > ${file_name}_tag.sam |& tee -a ${exec_log}
# var1=4 for the length of first part of tag (FLAG), for leading zero for correct sort
# var2 for the length of second tag (POS), for leading zero for correct sort


# sort according to BX:Z: TAG and remove duplicates
echo -e "\n\n################ Sort according to BX:Z: TAG\n\n" >> ${exec_log}
samtools sort -t "BX:Z:" -o ${file_name}_tag_sorted.bam ${file_name}_tag.sam # no need -h here
echo -e "\n\n################ Check the sorting\n\n" >> ${exec_log}
samtools view -h ${file_name}_tag_sorted.bam | head -30 >> ${exec_log}
echo -e "\n\n################ Remove duplicates\n\n" >> ${exec_log}
samtools view -h ${file_name}_tag_sorted.bam | awk '
BEGIN{COUNT=0}
/^@/{print $0 ; COUNT=int(COUNT+1) ; next}
{if(NR==COUNT){FLAG_0=$2 ; POS0_0=$4 ; POS1_0=int($4 + length($10) - 1) ; BEST=$5 ; LINE=$0 ; next}}
{FLAG_1=$2 ; POS0_1=$4 ; POS1_1=int($4 + length($10) - 1)}
{if(FLAG_0==FLAG_1 && POS0_0==POS0_1 && POS1_0==POS1_1){{if(BEST<$5){BEST=$5 ; LINE=$0}}}else{print LINE ; BEST=$5 ; LINE=$0}}
{FLAG_0=$2 ; POS0_0=$4 ; POS1_0=int($4 + length($10) - 1)}
END{print LINE}
' > ${file_name}_nodup.sam |& tee -a ${exec_log}
sed -i '/^$/d' ${file_name}_nodup.sam # remove empty lines. I do not know why I have an empty line after the header in the awk above
echo -e "\n\n################ Check the removal\n\n" >> ${exec_log}
cat ${file_name}_nodup.sam | head -30 >> ${exec_log}

# convertion into bam and sorting and indexing
echo -e "\n\n################ Convertion into bam and sorting and indexing\n\n" >> ${exec_log}
samtools view -bh -o ${file_name}_nodup.bam ${file_name}_nodup.sam |& tee -a ${exec_log}
samtools sort -o ${output_file} ${file_name}_nodup.bam |& tee -a ${exec_log}
samtools index ${output_file} |& tee -a ${exec_log} # create indexes

# check that no duplicates remains in TAG BX
echo -e "\n\n################ Check that no duplicates remains in TAG BX\n\n" >> ${exec_log}
echo -e "Number of lines with unique TAG:" >> ${exec_log}
awk '!/^@|^$/{print $NF}' ${file_name}_nodup.sam | sort | uniq -c | cut -d ' ' -f 1 | sort | uniq -c >> ${exec_log} # or rev file | cut -f 1
# !/^@|^$/ because one empty line after @ headers that I do see using samtools view I do not know why. Using Notepad+, I confirm the presence of the empty line
echo -e "\n\nTotal number of lines:" >> ${exec_log}
line_nb_after=$(samtools view ${output_file} | wc -l | cut -f1 -d' ')

# printing
LC_NUMERIC="en_US.UTF-8" # this is to have printf working for comma thousand separator
line_nb_before=$(samtools view ${bam} | wc -l | cut -f1 -d' ')
echo -e "\n\nNumber of sequences before removing duplicates: $(printf "%'d" ${line_nb_before})\n" >> ${log}
echo -e "Number of sequences after removing duplicates: $(printf "%'d" ${line_nb_after})\n" >> ${log}
echo -e "Ratio: " >> ${log}
echo -e $(printf '%.2f\n' $(echo " $line_nb_after / $line_nb_before " | bc -l)) >> ${log} # the number in '%.2f\n' is the number of decimals
echo -e "\n\n" >> ${log}
