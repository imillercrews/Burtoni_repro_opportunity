#### Reproduction synchronization burtoni

#### Downloading data from BaseSpace
##following directions from: https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-examples
## See also: https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
 
##enter directory:
#/stor/work/Hofmann/All_projects/A_burtoni_repro_opportunity/Synchronization/RNAseq/

## download basespace Linux
#run on lambcomp01
mkdir basespace
cd basespace
wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"

#give access to run  
chmod u+x ./bs

#authenticate basespace if needed
bs auth
#not needed already have auth file in /stor/home/imc/.basespace/default.cfg

#check that basespace is connected 
./bs whoami

#list projects
./bs list projects
#JA22396
#370823465 

### downloading samples
## save on stampede2
#create directory in work
cd ..
mkdir Data
cd Data

#burtoni sample ID: 370823465
#download project
../basespace/bs download project -i 370823465 -o basespace

## run fastqc
find ./basespace/ -name "*.fastq.gz" | xargs -n 1 fastqc -casava

## move files to fastqc folder
mkdir fastq
cd fastq
find ../basespace/ -name "*fastq*" | xargs -n 1 -I _ cp _ .

#remove basespace folder
cd ..
rm -r basespace

##run multiqc
multiqc ./fastq/

#### STAR
### reference
## Nile tilapia
## download reference
## use tilapia transcriptome
mkdir NileTilapia
cd NileTilapia
#annotation
wget http://ftp.ensembl.org/pub/release-105/gtf/oreochromis_niloticus/Oreochromis_niloticus.O_niloticus_UMD_NMBU.105.gtf.gz
#reference genome
wget http://ftp.ensembl.org/pub/release-105/fasta/oreochromis_niloticus/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa.gz

#gunzip
gunzip *.gz

#start new tmux
tmux new -s tilapia

# index reference
mkdir genomeOutput

STAR --runThreadN 6 --runMode genomeGenerate --genomeDir genomeOutput --genomeFastaFiles Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa  --sjdbGTFfile Oreochromis_niloticus.O_niloticus_UMD_NMBU.105.gtf --sjdbOverhang 99

#detach 
# Ctrl+b then d
#attach tmux a or tmux a -t tilapia

#type this in tmux and press 'enter' to get an email when done
#echo "tilapia done" | mail -s "Send an email with MAIL" imillercrews@utexas.edu
#rm -r 'tilapia done'

##Mapping .fastq file to reference
#need to list all fastq names
#cd Data
find ./fastq -name "*L001*.fastq.gz" > R1.txt

# concatenate fastq files
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_ME_R1_001.fastq.gz

done;

# move merged files to new folder
cd Data
mkdir merged_fastq
cd merged_fastq
mv ../fastq/*_ME_*.gz .

# gunzip files


## STAR alignment and count
mkdir bams
cd bams

STAR --genomeLoad LoadAndExit --genomeDir ../../NileTilapia/genomeOutput

for i in $(ls ../merged_fastq/*_ME_*.fastq | sed s/_[12].fq.gz// | sort -u)
do
    STAR --runMode alignReads --outSAMtype BAM Unsorted --genomeDir ../../NileTilapia/genomeOutput --readFilesIn ${i} --runThreadN 8 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2 --quantMode GeneCounts --outFileNamePrefix $(echo ${i} | while read F; do basename $F | rev | cut -c 22- | rev | sed s/$/-/; done)
done

STAR --genomeLoad Remove --genomeDir ../../NileTilapia/genomeOutput


echo "done" | mail -s "Send an email with MAIL" imillercrews@utexas.edu











############## trash ##############






## use burtoni new reference
#https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5321-6
#use novel annotation and novel reference transcripts
#"Time matters! Developmental shift in gene expression between the head and the trunk region of the cichlid fish Astatotilapia burtoni"
#anotation
cds
cd ref
mkdir burtoni
cd burtoni
#new names?
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-018-5321-6/MediaObjects/12864_2018_5321_MOESM5_ESM.txt
#new gtf
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-018-5321-6/MediaObjects/12864_2018_5321_MOESM3_ESM.gtf
#novel reference genome
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-018-5321-6/MediaObjects/12864_2018_5321_MOESM4_ESM.fasta
#reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/239/415/GCA_000239415.1_AstBur1.0/GCA_000239415.1_AstBur1.0_genomic.fna.gz
gunzip *.gz

#combine fasta files
mv 12864_2018_5321_MOESM4_ESM.fasta  12864_2018_5321_MOESM4_ESM.fna
cat *.fna > burtoni.concatenated.fasta








#test run
mkdir test
cd test
cp ../fastq/F2-1-J2-F_S170_ME_R1_001.fastq.gz .
cp ../fastq/H4-3-F3-F_S188_ME_R1_001.fastq.gz .
cp ../fastq/F2-1-J2-F_S170_L001_R1_001.fastq.gz .


STAR --runMode alignReads --outSAMtype BAM Unsorted --genomeDir ../NileTilapia/genomeOutput --outFileNamePrefix ./bams --readFilesIn `ls ./test/*_ME_*.fastq | sort | tr "\n" ","` --runThreadN 8 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2

$(ls -1 ./test/*_ME_*.fastq.gz | sort | tr "\n" ",")

STAR --runMode alignReads --outSAMtype BAM Unsorted --genomeDir ../NileTilapia/genomeOutput --readFilesIn ./test/F2-1-J2-F_S170_ME_R1_001.fastq --runThreadN 8 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2

STAR --runMode alignReads --outSAMtype BAM Unsorted --genomeDir ../NileTilapia/genomeOutput --readFilesIn ./test/F2-1-J2-F_S170_ME_R1_001.fastq --runThreadN 8 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2 --outFileNamePrefix SRR391535

#delete test files



STAR --runMode alignReads --outSAMtype BAM Unsorted --genomeDir ../NileTilapia/genomeOutput --readFilesIn $(ls ./test/*_ME_*.fastq | sort | tr "\n" ",") --runThreadN 8 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2 --outFileNamePrefix $(find ./test/ -type f -name "*.fastq" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort)




mkdir bams
cd bams

STAR --genomeLoad LoadAndExit --genomeDir ../../NileTilapia/genomeOutput

for i in $(ls ../test/*_ME_*.fastq | sed s/_[12].fq.gz// | sort -u)
do
    STAR --runMode alignReads --outSAMtype BAM Unsorted --genomeDir ../../NileTilapia/genomeOutput --readFilesIn ${i} --runThreadN 8 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2 --outFileNamePrefix $(echo ${i} | while read F; do basename $F | rev | cut -c 22- | rev | sed s/$/-/; done)
done

STAR --genomeLoad Remove --genomeDir ../../NileTilapia/genomeOutput







for i in $(ls ../test/*_ME_*.fastq | sed s/_[12].fq.gz// | sort -u)
do
    echo ${i} | while read F; do basename $F | rev | cut -c 22- | rev | sed s/$/-/; done
done





STAR --genomeLoad LoadAndExit --genomeDir index.150

for i in $(ls raw_data | sed s/_[12].fq.gz// | sort -u)
do
    STAR [...]
done

STAR --genomeLoad Remove --genomeDir index.150




for i in $(ls 1_raw_data/*_R1_* | sort -u); do echo STAR --genomeDir /home/pahib/RNA_SEQ_Pipeline/Reference_genome/Drosophila_STAR/ \
--readFilesIn 1_raw_data/${i} 1_raw_data/${i/_R1_/_R2_} \
--runThreadN 20 --outFileNamePrefix 3_aligned/${i/_R1_001.fastq.gz/} \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--sjdbGTFfile /home/pahib/RNA_SEQ_Pipeline/Reference_genome/Drosophila_gtf/Drosophila_melanogaster.BDGP6.32.104.chr.gtf \
--readFilesCommand gunzip -c ; done;











#run on all merged files
STAR --runMode alignReads --outSAMtype BAM Unsorted --genomeDir ./NileTilapia --outFileNamePrefix ./bams --readFilesIn `ls ./fastq/*_ME_*.fastq.gz | sort | tr "\n" ","` --runThreadN 8 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2 --readFilesCommand gunzip -c






### Count data
cd Data
mkdir count
cd count

featureCounts -t exon -g gene_id -a /path/to/annotation.gtf -o /path/to/output.txt /path/to/mapping_results_SE.bam

featureCounts -t exon -g gene_id -a ../../NileTilapia/Oreochromis_niloticus.O_niloticus_UMD_NMBU.105.gtf -o burtoni_repro_counts.txt ../bams/*.bam


featureCounts -T 4 -s 2 \
  -a /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/Homo_sapiens.GRCh38.92.gtf \
  -o ~/unix_lesson/rnaseq/results/counts/Mov10_featurecounts.txt \
  ~/unix_lesson/rnaseq/results/STAR/bams/*.out.bam










STAR --genomeLoad LoadAndExit --genomeDir ../../NileTilapia/genomeOutput

for i in $(ls ./*_ME_*.fastq | sed s/_[12].fq.gz// | sort -u)
do
    STAR --runMode alignReads --outSAMtype BAM Unsorted --genomeDir ../../NileTilapia/genomeOutput --readFilesIn ${i} --runThreadN 8 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 2 --quantMode GeneCounts --outFileNamePrefix $(echo ${i} | while read F; do basename $F | rev | cut -c 22- | rev | sed s/$/-/; done)
done

STAR --genomeLoad Remove --genomeDir ../../NileTilapia/genomeOutput





### mark duplicates
A1-1-A3-M-Aligned.out.bam 

java -jar picard.jar MarkDuplicates --ASSUME_SORT_ORDER null \
      I= ../bams/A1-1-A3-M-Aligned.out.bam  \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt 

java -jar picard.jar MarkDuplicates --ASSUME_SORT_ORDER null \
      I= ../bams/A1-1-A3-M-Aligned.out.bam  \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt \
	ASSUME_SORTED=true


java -jar $PICARD MarkDuplicates \
      I= ../bams/A1-1-A3-M-Aligned.out.bam  \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt





