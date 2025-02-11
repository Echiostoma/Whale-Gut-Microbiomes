#!/bin/bash
#SBATCH --cpus-per-task=16 --mem-per-cpu=16g
mamba activate bowtie2
bowtie2 -x whale -1 Hi1A_L1_L2_R1_val_1.fq.gz -2 Hi1A_L1_L2_R2_val_2.fq.gz -S Hi1Amapped.sam --very-sensitive-local --un-conc-gz Hi1A.fastq.gz

bowtie2 -x whale -1 Hi1D_L1_L2_R1_val_1.fq.gz -2 Hi1D_L1_L2_R2_val_2.fq.gz -S Hi1Dmapped.sam --very-sensitive-local --un-conc-gz Hi1D.fastq.gz

bowtie2 -x whale -1 MN1_L1_L2_R1_val_1.fq.gz -2 MN1_L1_L2_R2_val_2.fq.gz -S MN1mapped.sam --very-sensitive-local --un-conc-gz MN1.fastq.gz

bowtie2 -x whale -1 MN3_L1_L2_R1_val_1.fq.gz -2 MN3_L1_L2_R2_val_2.fq.gz -S MN3mapped.sam --very-sensitive-local --un-conc-gz MN3.fastq.gz

bowtie2 -x whale -1 MN4_L1_L2_R1_val_1.fq.gz -2 MN4_L1_L2_R2_val_2.fq.gz -S MN4mapped.sam --very-sensitive-local --un-conc-gz MN4.fastq.gz


#!/bin/bash
#SBATCH --cpus-per-task=16 --mem-per-cpu=24g
mamba activate kaiju
kaiju-makedb -s nr

#!/bin/bash
#SBATCH --cpus-per-task=16 --mem-per-cpu=16g
#SBATCH -p sapphire # Partition to submit to
#SBATCH -c 16 # Number of threads
mamba activate kaiju
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK kaiju-multi -z 16 -t nodes.dmp -f /n/holyscratch01/girguis_lab/Lab/rbos/kaijudb/nr/kaiju_db_nr.fmi -i Hi1A.fastq.1.gz,Hi1d.1.fastq,MN1.1.fastq,MN3.1.fastq,MN4.1.fastq  -j Hi1A.fastq.2.gz,Hi1d.2.fastq,MN1.2.fastq,MN3.2.fastq,MN4.2.fastq -o Hi1A.out,Hi1D.out,MN1.out,MN3.out,MN4.out -E .0001

kaiju2table -t nodes.dmp -n names.dmp -l phylum,family -m 0.5 -o kaiju_summary.tsv Hi1A.out Hi1D.out MN1.out MN3.out MN4.out

#!/bin/bash
#SBATCH -c 64 # Number of threads
#SBATCH -t 1-0:00:00
#SBATCH -p sapphire # Partition to submit to
#SBATCH --mem-per-cpu=8g #Memory per cpu
#SBATCH --cpus-per-task=64
mamba activate spades
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK metaspades.py -1 Hi1A.1.fastq.gz -2 Hi1A.2.fastq.gz -t 64 -m 512 -o Hi1A_Assembly

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK metaspades.py -1 Hi1d.1.fastq.gz -2 Hi1d.2.fastq.gz -t 64 -m 512 -o Hi1D_Assembly

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK metaspades.py -1 MN1.1.fastq.gz -2 MN1.2.fastq.gz -t 64 -m 512 -o MN1_Assembly

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK metaspades.py -1 MN3.1.fastq.gz -2 MN3.2.fastq.gz -t 64 -m 512 -o MN3_Assembly

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK metaspades.py -1 MN4.1.fastq.gz -2 MN4.2.fastq.gz -t 64 -m 512 -o MN4_Assembly

mamba activate bowtie2
bowtie2 -x Hi1A -1 Hi1A.1.fastq.gz -2 Hi1A.2.fastq.gz -S Hi1A_Assembly_mapped.sam --very-sensitive-local

bowtie2 -x Hi1D -1 Hi1D.1.fastq.gz -2 Hi1D.2.fastq.gz -S Hi1D_Assembly_mapped.sam --very-sensitive-local

bowtie2 -x MN1 -1 MN1.1.fastq.gz -2 MN1.2.fastq.gz -S MN1_Assembly_mapped.sam --very-sensitive-local

bowtie2 -x MN3 -1 MN3.1.fastq.gz -2 MN3.2.fastq.gz -S MN3_Assembly_mapped.sam --very-sensitive-local

bowtie2 -x MN4 -1 MN4.1.fastq.gz -2 MN4.2.fastq.gz -S MN4_Assembly_mapped.sam --very-sensitive-local

#!/bin/bash
#SBATCH -c 16 # Number of threads
#SBATCH -t 1-0:00:00
#SBATCH -p sapphire # Partition to submit to
#SBATCH --mem-per-cpu=8g #Memory per cpu
#SBATCH --cpus-per-task=16
mamba activate samtools
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK samtools sort Hi1A_Assembly_mapped.sam -o Hi1A_Assembly_mapped_sorted.bam --threads 16

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK samtools sort Hi1D_Assembly_mapped.sam -o Hi1D_Assembly_mapped_sorted.bam --threads 16

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK samtools sort MN1_Assembly_mapped.sam -o MN1_Assembly_mapped_sorted.bam --threads 16

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK samtools sort MN3_Assembly_mapped.sam -o MN3_Assembly_mapped_sorted.bam --threads 16

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK samtools sort MN4_Assembly_mapped.sam -o MN4_Assembly_mapped_sorted.bam --threads 16

mamba activate metabat2
runMetaBat.sh Hi1A_Assembly.fasta_300_.fasta Hi1A_Assembly_mapped_sorted.bam

runMetaBat.sh Hi1D_Assembly.fasta_300_.fasta Hi1D_Assembly_mapped_sorted.bam

runMetaBat.sh MN1_Assembly.fasta_300_.fasta MN1_Assembly_mapped_sorted.bam

runMetaBat.sh MN3_Assembly.fasta_300_.fasta MN3_Assembly_mapped_sorted.bam

runMetaBat.sh MN4_Assembly.fasta_300_.fasta MN4_Assembly_mapped_sorted.bam

metabat2 -i Hi1A_Assembly.fasta_300_.fasta -a Hi1A_depth.txt -o bins_dir/bin

metabat2 -i Hi1D_Assembly.fasta_300_.fasta -a Hi1D_depth.txt -o bins_dir/bin

metabat2 -i MN1_Assembly.fasta_300_.fasta -a MN1_depth.txt -o bins_dir/bin

metabat2 -i MN3_Assembly.fasta_300_.fasta -a MN3_depth.txt -o bins_dir/bin

metabat2 -i MN4_Assembly.fasta_300_.fasta -a MN4_depth.txt -o bins_dir/bin

#!/bin/bash
#SBATCH -c 16 # Number of threads
#SBATCH -t 1-0:00:00
#SBATCH -p sapphire # Partition to submit to
#SBATCH --mem-per-cpu=8g #Memory per cpu
#SBATCH --cpus-per-task=16
mamba activate anvio-8
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK for f in *.fasta; do anvi-gen-contigs-database -T 16 -f $f -n ${f} -o ${f}_out.db; done

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK for f in *.db; do anvi-run-ncbi-cogs -c $f -T 16; done

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK for f in *.db; do anvi-run-kegg-kofams -c $f -T 16; done

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK for f in *.db; do anvi-run-pfams -c $f -T 16; done

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK for f in *.db; do anvi-run-cazymes -c $f -T 16; done

#!/bin/bash
#SBATCH --cpus-per-task=16 --mem-per-cpu=4g
mamba activate anvio-8
###for all files (f) with the .db suffix (*.db) in the working directory, export (anvi-export-functions) annotation sources (COG20_CATEGORY, COG20_FUNCTION, KOfam, Pfam) for each contig database (-c) with the pasted contig database name ($f) to a new .txt file (-o) with the pasted ($f) contig database name###
for f in *.db; do anvi-export-functions -c $f -o ${f}_COG20Category.txt --annotation-sources COG20_CATEGORY; done
for f in *.db; do anvi-export-functions -c $f -o ${f}_COG20Function.txt --annotation-sources COG20_FUNCTION; done
for f in *.db; do anvi-export-functions -c $f -o ${f}_KOfam.txt --annotation-sources KOfam; done
for f in *.db; do anvi-export-functions -c $f -o ${f}_Pfam.txt --annotation-sources Pfam; done
for f in *.db; do anvi-export-functions -c $f -o ${f}_Cazyme.txt --annotation-sources CAZyme; done

###for all files (f) with the .db suffix (*.db) in the working directory, export prodigal (--gene-caller prodigal) designated gene calls (anvi-export-gene-calls) for each contig database (-c) to a new text file (-o) with the pasted contig database name###
for f in *.db; do anvi-export-gene-calls -c $f --gene-caller prodigal -o ${f}_AllGeneCalls.txt; done
for f in *.db; do anvi-get-sequences-for-gene-calls -c $f -o ${f}_SequencesForGeneCalls.txt; done

###for all files (f) with the COG20Function.txt suffix, run the awk command to filter column 5 ($5) to include genes with an e-value <= 0.0001 and export (>) to a new file named ${f}_COG20FunctionEvalue_filtered.txt###
for f in *COG20Function.txt; do awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' $f >${f}_COG20FunctionEvalue_filtered.txt;done

###for all files (f) with the KOfam.txt suffix, run the awk command to filter column 5 ($5) to include genes with an e-value <= 0.0001 and export (>) to a new file named ${f}_KOfamEvalue_filtered.txt###
for f in *KOfam.txt; do awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' $f >${f}_KOfamEvalue_filtered.txt;done

###for all files (f) with the Pfam.txt suffix, run the awk command to filter column 5 ($5) to include genes with an e-value <= 0.0001 and export (>) to a new file named ${f}_PfamEvalue_filtered.txt###
for f in *Pfam.txt; do awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' $f >${f}_PfamEvalue_filtered.txt;done

###for all files (f) with the *KOfamEvalue_filtered.txt suffix, cut all genes with e-values <=0.0001 from the previous step, sort them alphabetically (sort), retain only uniq values (uniq -c) while counting the number of appearances of each gene and export (>) to a new file named ${f}KOfamCounts.txt###
for f in *KOfamEvalue_filtered.txt; do cut -f4 $f| sort | uniq -c > ${f}KOfamCounts.txt; done

###for all files (f) with the *PfamEvalue_filtered.txt suffix, cut all genes with e-values <=0.0001 from the previous step, sort them alphabetically (sort), retain only uniq values (uniq -c) while counting the number of appearances of each gene and export (>) to a new file named ${f}PfamCounts.txt###
for f in *PfamEvalue_filtered.txt; do cut -f4 $f| sort | uniq -c > ${f}PfamCounts.txt; done

###for all files (f) with the *COG20FunctionEvalue_filtered.txt suffix, cut all genes with e-values <=0.0001 from the previous step, sort them alphabetically (sort), retain only uniq values (uniq -c) while counting the number of appearances of each gene and export (>) to a new file named ${f}COG20FunctionCounts.txt###
for f in *COG20FunctionEvalue_filtered.txt; do cut -f4 $f| sort | uniq -c > ${f}COG20FunctionCounts.txt; done


###for all files (f) with the suffix KOfamCounts.txt; run the sed command to remove white space introduced by awk and export (>) to a new file ${f}KOfamCountsReformatted.txt###
for f in *KOfamCounts.txt; do sed 's/^ *//g' $f > ${f}KOfamCountsReformatted.txt;done

###for all files (f) with the suffix KOfamCountsReformatted.txt; run the sed command to place a semicolon between gene name and counts for that gene and export (>) to a new file ${f}KOfamCountsReformattedTwo.txt###
for f in *KOfamCountsReformatted.txt; do sed 's/ /;/' $f > ${f}KOfamCountsReformattedTwo.txt; done

###for all files (f) with the suffix KOfamCountsReformattedTwo.txt, do awk to split columns one and two into two new cells at the semi colon by gene name ($1 column1) and gene counts ($2 column2) to a new file called ${f}KOfamFinalCounts.tsv###
for f in *KOfamCountsReformattedTwo.txt; do awk -F ";" '{n=split($1,a,";");for (i=1;i<=n;i++) print $2"\t"a[i]}' $f > ${f}KOfamFinalCounts.tsv; done

###for all files (f) with the suffix PfamCounts.txt; run the sed command to remove white space introduced by awk and export (>) to a new file ${f}PfamCountsReformatted.txt###
for f in *PfamCounts.txt; do sed 's/^ *//g' $f > ${f}PfamCountsReformatted.txt;done

###for all files (f) with the suffix PfamCountsReformatted.txt; run the sed command to place a semicolon between gene name and counts for that gene and export (>) to a new file ${f}PfamCountsReformattedTwo.txt###
for f in *PfamCountsReformatted.txt; do sed 's/ /;/' $f > ${f}PfamCountsReformattedTwo.txt; done

###for all files (f) with the suffix PfamCountsReformattedTwo.txt, do awk to split columns one and two into two new cells at the semi colon by gene name ($1 column1) and gene counts ($2 column2) to a new file called ${f}PfamFinalCounts.tsv###
for f in *PfamCountsReformattedTwo.txt; do awk -F ";" '{n=split($1,a,";");for (i=1;i<=n;i++) print $2"\t"a[i]}' $f > ${f}PfamFinalCounts.tsv; done

###for all files (f) with the suffix COG20FunctionCounts.txt; run the sed command to remove white space introduced by awk and export (>) to a new file ${f}COG20FunctionCountsReformatted.txt###
for f in *COG20FunctionCounts.txt; do sed 's/^ *//g' $f > ${f}COG20FunctionCountsReformatted.txt;done

###for all files (f) with the suffix COG20FunctionCountsReformatted.txt; run the sed command to place a semicolon between gene name and counts for that gene and export (>) to a new file ${f}COG20FunctionCountsReformattedTwo.txt###
for f in *COG20FunctionCountsReformatted.txt; do sed 's/ /;/' $f > ${f}COG20FunctionCountsReformattedTwo.txt; done

###for all files (f) with the suffix COG20FunctionCountsReformattedTwo.txt, do awk to split columns one and two into two new cells at the semi colon by gene name ($1 column1) and gene counts ($2 column2) to a new file called ${f}COG20FunctionFinalCounts.tsv###
for f in *COG20FunctionCountsReformattedTwo.txt; do awk -F ";" '{n=split($1,a,";");for (i=1;i<=n;i++) print $2"\t"a[i]}' $f > ${f}COG20FunctionFinalCounts.tsv; done


###list each gene file with suffix *KOfamFinalCounts.tsv and cut the gene names from column 1 of each file and export to new file called all_genes.txt
for gene_file in `ls *KOfamFinalCounts.tsv`; do
  cut -f1 ${gene_file} >> all_genes.txt
done


###sort all of the genes alphabetically, and only retains unique occurences of each gene and export to new file called uniq_KOfam.txt###
sort -u all_genes.txt > uniq_KOfam.txt


###the next portion of this script is run a total of three times, once for each annotation source, so only this first occurrence is annotated for KOfams here###
###set the environment we are working in###
set -e


###establish the variable called uniq_genes_file and define it as uniq_KOfam.txt###
uniq_genes_file="uniq_KOfam.txt"


###establish the variable output_file and define it as combined-outputKOfams.tsv###
output_file="combined-outputKOfams.tsv"

###append/create the beginning of the output file with column name added, which contains all unique gene names that we identified above in the first column###
cat <( printf "function\n" ) ${uniq_genes_file} > ${output_file}


###set to break only on a new line only###
IFS=$'\n'

###loop through all of our input files with the suffix KOfamFinalCounts.tsv###
for curr_file in $(ls *KOfamFinalCounts.tsv); do

###create a sample name without the filename extension###

curr_sample_name=$(echo ${curr_file} | sed 's/.tsv//')

###append the current sample name to a temporary file that we build over time###
printf "${curr_sample_name}\n" > building-col.tmp

###loop through unique gene names in the current sample file to append counts or '0' if that gene name is not matched exactly with the grep command below###
for curr_gene in $(cat ${uniq_genes_file}); do


if grep -m 1 -F -q "${curr_gene}" ${curr_file}; then

grep -m 1 -F "${curr_gene}" ${curr_file} | cut -f 2 >> building-col.tmp

else

  printf "0\n" >> building-col.tmp

fi

done

###add to our output file that is being built over time and remove the current sample file until the last samile file is reached and processed###
paste ${output_file} building-col.tmp > building-output.tmp
mv building-output.tmp ${output_file}
rm building-col.tmp

done

###delete all files with the following suffixes###
find . -name "*KOfamCounts.txt" -type f -delete
find . -name "*KOfamEvalue_filtered.txt" -type f -delete
find . -name "*Reformatted.txt" -type f -delete
find . -name "*ReformattedTwo.txt" -type f -delete
find . -name "*KOfamFinalCounts.tsv" -type f -delete




###to create gene catalogs with pfam annotations###
for gene_file in `ls *PfamFinalCounts.tsv`; do
  cut -f1 ${gene_file} >> all_Pfam.txt
done

sort -u all_Pfam.txt > uniq_Pfams.txt


set -e

uniq_genes_file="uniq_Pfams.txt"

output_file="combined-outputPfams.tsv"


cat <( printf "function\n" ) ${uniq_genes_file} > ${output_file}


IFS=$'\n'


for curr_file in $(ls *PfamFinalCounts.tsv); do


curr_sample_name=$(echo ${curr_file} | sed 's/.tsv//')


printf "${curr_sample_name}\n" > building-col.tmp


for curr_gene in $(cat ${uniq_genes_file}); do


if grep -m 1 -F -q "${curr_gene}" ${curr_file}; then

grep -m 1 -F "${curr_gene}" ${curr_file} | cut -f 2 >> building-col.tmp

else

  printf "0\n" >> building-col.tmp

fi

done


paste ${output_file} building-col.tmp > building-output.tmp
mv building-output.tmp ${output_file}
rm building-col.tmp

done

###delete all files with the following suffixes###
find . -name "*PfamCounts.txt" -type f -delete
find . -name "*PfamEvalue_filtered.txt" -type f -delete
find . -name "*Reformatted.txt" -type f -delete
find . -name "*ReformattedTwo.txt" -type f -delete
find . -name "*PfamFinalCounts.tsv" -type f -delete



###to create gene catalogs with COG20 Functions###
for gene_file in `ls *COG20FunctionFinalCounts.tsv`; do
  cut -f1 ${gene_file} >> all_COG20Function.txt
done

sort -u all_COG20Function.txt > uniq_COG20Function.txt


set -e

uniq_genes_file="uniq_COG20Function.txt"

output_file="combined-outputCOG20Function.tsv"


cat <( printf "function\n" ) ${uniq_genes_file} > ${output_file}


IFS=$'\n'


for curr_file in $(ls *COG20FunctionFinalCounts.tsv); do


curr_sample_name=$(echo ${curr_file} | sed 's/.tsv//')


printf "${curr_sample_name}\n" > building-col.tmp


for curr_gene in $(cat ${uniq_genes_file}); do


if grep -m 1 -F -q "${curr_gene}" ${curr_file}; then

grep -m 1 -F "${curr_gene}" ${curr_file} | cut -f 2 >> building-col.tmp

else

  printf "0\n" >> building-col.tmp

fi

done


paste ${output_file} building-col.tmp > building-output.tmp
mv building-output.tmp ${output_file}
rm building-col.tmp

done

###delete all files with the following suffixes###
find . -name "*COG20FunctionCounts.txt" -type f -delete
find . -name "*COG20FunctionEvalue_filtered.txt" -type f -delete
find . -name "*Reformatted.txt" -type f -delete
find . -name "*ReformattedTwo.txt" -type f -delete
find . -name "*COG20FunctionFinalCounts.tsv" -type f -delete

#!/bin/bash
#SBATCH -t 1-0:00:00
#SBATCH -p sapphire # Partition to submit to
#SBATCH --mem-per-cpu=8g #Memory per cpu
#SBATCH --cpus-per-task=16
mamba activate checkm2
checkm2 database --download --path /n/holyscratch01/girguis_lab/ryan_bos/whales

#!/bin/bash
#SBATCH -c 16 # Number of threads
#SBATCH -t 1-0:00:00
#SBATCH -p sapphire # Partition to submit to
#SBATCH --mem-per-cpu=8g #Memory per cpu
#SBATCH --cpus-per-task=16
mamba activate checkm2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK checkm2 predict --threads 16 -x .fa --input /n/holyscratch01/girguis_lab/Lab/rbos/Hi1ANoCheckMbins_dir --output-directory /n/holyscratch01/girguis_lab/Lab/rbos/Hi1ANoCheckMbins_dir/CheckM2results --database_path /n/holyscratch01/girguis_lab/ryan_bos/whales/CheckM2_database/uniref100.KO.1.dmnd
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK checkm2 predict --threads 16 -x .fa --input /n/holyscratch01/girguis_lab/Lab/rbos/Hi1DNoCheckMbins_dir --output-directory /n/holyscratch01/girguis_lab/Lab/rbos/Hi1DNoCheckMbins_dir/CheckM2results --database_path /n/holyscratch01/girguis_lab/ryan_bos/whales/CheckM2_database/uniref100.KO.1.dmnd
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK checkm2 predict --threads 16 -x .fa --input /n/holyscratch01/girguis_lab/Lab/rbos/MN1NoCheckMbins_dir --output-directory /n/holyscratch01/girguis_lab/Lab/rbos/MN1NoCheckMbins_dir/CheckM2results --database_path /n/holyscratch01/girguis_lab/ryan_bos/whales/CheckM2_database/uniref100.KO.1.dmnd
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK checkm2 predict --threads 16 -x .fa --input /n/holyscratch01/girguis_lab/Lab/rbos/MN3NoCheckMbins_dir --output-directory /n/holyscratch01/girguis_lab/Lab/rbos/MN3NoCheckMbins_dir/CheckM2results --database_path /n/holyscratch01/girguis_lab/ryan_bos/whales/CheckM2_database/uniref100.KO.1.dmnd
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK checkm2 predict --threads 16 -x .fa --input /n/holyscratch01/girguis_lab/Lab/rbos/MN4NoCheckMbins_dir --output-directory /n/holyscratch01/girguis_lab/Lab/rbos/MN4NoCheckMbins_dir/CheckM2results --database_path /n/holyscratch01/girguis_lab/ryan_bos/whales/CheckM2_database/uniref100.KO.1.dmnd



#!/bin/bash
#SBATCH -c 32 # Number of threads
#SBATCH -p sapphire # Partition to submit to
#SBATCH --mem-per-cpu=8g #Memory per cpu
#SBATCH --cpus-per-task=32
mamba activate gtdbtk-2.4.0
export GTDBTK_DATA_PATH="/n/holyscratch01/girguis_lab/Lab/rbos/assemblies/release220"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK gtdbtk classify_wf --extension .fa --genome_dir /n/holyscratch01/girguis_lab/Lab/rbos/MN1NoCheckMbins_dir --out_dir /n/holyscratch01/girguis_lab/Lab/rbos/MN1NoCheckMbins_dir/GTDBtkResults --cpus 32 --mash_db /n/holyscratch01/girguis_lab/Lab/rbos/assemblies/gtdb-tk-r220.msh

#!/bin/bash
#SBATCH -c 32 # Number of threads
#SBATCH -t 1-0:00:00
#SBATCH -p sapphire # Partition to submit to
#SBATCH --mem-per-cpu=4g #Memory per cpu
#SBATCH --cpus-per-task=32
mamba activate gtotree
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -c $SLURM_CPUS_PER_TASK GToTree -a bacillota.txt -f mags.txt -H Firmicutes -j 32 -o Firmicutes
srun -c $SLURM_CPUS_PER_TASK GToTree -a bacteroides.txt -f mags.txt -H Bacteroidetes -j 32 -o Bacteroidetes
srun -c $SLURM_CPUS_PER_TASK GToTree -a actinobacteria.txt -f mags.txt -H Actinobacteria -j 32 -o Actinobacteria
srun -c $SLURM_CPUS_PER_TASK GToTree -a VerrucomicrobiotaANDspirochaetota.txt -f mags.txt -H Bacteria -j 32 -o VerruocomicrobiotaANDspirochaetota

mamba activate smash
gunzip genbank-2022.03-bacteria-k31.zip
ls -lh genbank-2022.03-bacteria-k31.zip
gunzip genbank-2022.03-bacteria.lineages.csv.gz
sourmash sketch dna -p k=31,abund Hi1A*.gz -o Hi1A.sig.gz --name Hi1A
sourmash sketch dna -p k=31,abund Hi1d*.gz -o Hi1D.sig.gz --name Hi1D
sourmash sketch dna -p k=31,abund MN1*.gz -o MN1.sig.gz --name MN1
sourmash sketch dna -p k=31,abund MN3*.gz -o MN3.sig.gz --name MN3
sourmash sketch dna -p k=31,abund MN4*.gz -o MN4.sig.gz --name MN4

sourmash sketch dna Hi1A*.fa --name-from-first
sourmash sketch dna Hi1D*.fa --name-from-first
sourmash sketch dna MN1*.fa --name-from-first
sourmash sketch dna MN3*.fa --name-from-first
sourmash sketch dna MN4*.fa --name-from-first

sourmash gather Hi1A.sig.gz genbank-2022.03-bacteria-k31.zip --save-matches Hi1AGenbankMatches.zip
sourmash gather Hi1A.sig.gz Hi1A*.sig Hi1AGenbankMatches.zip -o Hi1AGenbankandMAGs.csv

sourmash gather Hi1D.sig.gz genbank-2022.03-bacteria-k31.zip --save-matches Hi1DGenbankMatches.zip
sourmash gather Hi1D.sig.gz Hi1D*.sig Hi1DGenbankMatches.zip -o Hi1DGenbankandMAGs.csv

sourmash gather MN1.sig.gz genbank-2022.03-bacteria-k31.zip --save-matches MN1GenbankMatches.zip
sourmash gather MN1.sig.gz MN1*.sig MN1GenbankMatches.zip -o MN1GenbankandMAGs.csv

sourmash gather MN3.sig.gz genbank-2022.03-bacteria-k31.zip --save-matches MN3GenbankMatches.zip
sourmash gather MN3.sig.gz MN3*.sig MN3GenbankMatches.zip -o MN3GenbankandMAGs.csv

sourmash gather MN4.sig.gz genbank-2022.03-bacteria-k31.zip --save-matches MN4GenbankMatches.zip
sourmash gather MN4.sig.gz MN4*.sig MN4GenbankMatches.zip -o MN4GenbankandMAGs.csv

mamba activate qiime2
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Manifest.txt --output-path Whales.qza --input-format PairedEndFastqManifestPhred33V2
qiime demux summarize --i-data Whales.qza --o-visualization Whales.qzv
qiime vsearch merge-pairs --i-demultiplexed-seqs Whales.qza --o-merged-sequences Whales.joined.qza --o-unmerged-sequences Whales --verbose
qiime demux summarize --i-data Whales.joined.qza --o-visualization Whales.joined.qzv --verbose
qiime quality-filter q-score \
    --i-demux Whales.joined.qza \
    --o-filtered-sequences Whales-demux-filtered.qza \
    --o-filter-stats Whales-demux-filter-stats.qza

qiime demux summarize --i-data Whales-demux-filtered.qza --o-visualization Whales-demux-filtered.qzv

qiime deblur denoise-16S --i-demultiplexed-seqs Whales-demux-filtered.qza --p-trim-length 225 --p-sample-stats --p-min-reads 0 --p-jobs-to-start 4 --o-table Whaleoutputtable --o-representative-sequences Whalesrepresentativesequences.qza --o-stats Whaleoutputstats
qiime feature-table tabulate-seqs --i-data Whalesrepresentativesequences.qza --o-visualization Whalesrepresentativesequences.qzv
qiime deblur visualize-stats --i-deblur-stats Whaleoutputstats.qza --o-visualization Whaleoutputstats.qzv
qiime feature-classifier classify-sklearn --i-classifier silva-138.1-ssu-nr99-515f-806r-classifier.qza --i-reads Whalesrepresentativesequences.qza --o-classification Whalereptaxonomy --p-reads-per-batch 10000 --verbose
qiime metadata tabulate --m-input-file Whalereptaxonomy.qza --o-visualization whalereptaxonomy.qzv
qiime taxa barplot --i-table Whalereptaxonomy.qza --i-taxonomy Whalereptaxonomy --o-visualization Whalerepfinal.qzv --m-metadata-file Manifest.qza
qiime taxa collapse \
  --i-table Whaleoutputtable.qza \
  --i-taxonomy Whalereptaxonomy.qza \
  --p-level 7 \
  --o-collapsed-table phyla-table.qza

qiime tools export --input-path Whalereptaxonomy.qza --output-path /Users/girguislab/kmers/qiimewhales/
biom add-metadata -i feature-table.biom -o table-with-taxonomy.biom --observation-metadata-fp taxonomy.tsv --sc-separated taxonomy
