source activate qiime2-2018.11

###################################################
###################################################
####### Demultiplexing by target gene #############
###################################################
###################################################

# Directories
DATADIR="fastq_files"
Primers="/home/olivierl/apps/Metagenetics_Projects/Primers/Primers_16S_18S.txt"
mkdir Demultiplexed_fastq_files
Output="Demultiplexed_fastq_files"
mkdir $Output/fastq_16S
mkdir $Output/fastq_18S

# Demultiplexing
for i in $DATADIR/*_R1_001.fastq.gz
do
R1=$i;
#replace "_R1_" with "_R2_" to get R2 file name
R2=${R1%_R1_001.fastq.gz}_R2_001.fastq.gz;
SAMPLE_F=`basename ${R1};`
SAMPLE_R=`basename ${R2};`
while read -r Primer_name Forward_primers Reverse_primers; do
cutadapt -g X$Forward_primers -G X$Reverse_primers -o $Output/$Primer_name$SAMPLE_F -p $Output/$Primer_name$SAMPLE_R $R1 $R2 -j 0 --no-indels --discard-untrimmed --overlap 17
done<$Primers
done > demultiplexing_summary.txt

# Moving demultiplexed fastq files into their relevant folders
mv $Output/16S* $Output/fastq_16S
mv $Output/18S* $Output/fastq_18S


###################################################
###################################################
################ 16S data #########################
###################################################
###################################################

# Directories
Database_fasta='/home/olivierl/apps/Metagenetics_Projects/Databases/SILVA_132_QIIME_release/silva_132_99_16S.fna'
Database_txt='/home/olivierl/apps/Metagenetics_Projects/Databases/SILVA_132_QIIME_release/16S_taxonomy_7_levels.txt'
Metadata_file='/home/olivierl/lus/fish_farms/Metadata_16S.txt'
mkdir /home/olivierl/lus/fish_farms/Results_16S_Dada2
Outputs='/home/olivierl/lus/fish_farms/Results_16S_Dada2'
fastqManifest='/home/olivierl/lus/fish_farms/Fastq_manifest_16S.csv'

#Importing files into qiime format
printf "\nImporting fastq files to qiime2\n" 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $fastqManifest \
  --output-path $Outputs/paired_end_demux.qza \
  --source-format PairedEndFastqManifestPhred33

#Visualization of sequence quality
qiime demux summarize \
  --i-data $Outputs/paired_end_demux.qza \
  --o-visualization $Outputs/paired_end_demux.qzv \
  --verbose 
  
# Quality filtering, denoising and chimera filtering using DADA2
printf "\nStarting DADA2 process\n" 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $Outputs/paired_end_demux.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 228 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 216 \
  --p-chimera-method consensus \
  --p-min-fold-parent-over-abundance 1 \
  --p-n-threads 0 \
  --o-representative-sequences $Outputs/rep_seqs_dada2.qza \
  --o-table $Outputs/table_dada2.qza \
  --o-denoising-stats $Outputs/stats-dada2.qza \
  --verbose


printf "\nCreating visualisation files of DADA2 outputs\n" 
# Looking at the feature table and feature data summary
qiime feature-table summarize \
  --i-table $Outputs/table_dada2.qza \
  --o-visualization $Outputs/table_dada2.qzv \
  --m-sample-metadata-file $Metadata_file

qiime feature-table tabulate-seqs \
  --i-data $Outputs/rep_seqs_dada2.qza \
  --o-visualization $Outputs/rep_seqs_dada2.qzv
  
qiime metadata tabulate \
   --m-input-file $Outputs/stats-dada2.qza \
   --o-visualization $Outputs/stats-dada2.qzv
   
   
# Training a Naive Bayes classifier - assigning taxonomy from different databases
# Importing the fasta and taxonomy files
printf "\nInitiating taxonomic assignment\n" 
printf "\nImporting reference database\n" 

qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path $Database_fasta \
   --output-path $Outputs/database_seq.qza

 qiime tools import \
   --type 'FeatureData[Taxonomy]' \
   --source-format HeaderlessTSVTaxonomyFormat \
   --input-path $Database_txt \
   --output-path $Outputs/database_taxonomy.qza
   
# Extracting the sequence portion we have targeted
printf "\nExtracting the sequence portion we have targeted\n" 

 qiime feature-classifier extract-reads \
   --i-sequences $Outputs/database_seq.qza \
   --p-f-primer CCTACGGGNGGCWGCAG \
   --p-r-primer GACTACHVGGGTSTCTAATCC \
   --o-reads $Outputs/database_seq.qza
   
# Trainning the classifier
printf "\nTrainning the classifier\n" 

 qiime feature-classifier fit-classifier-naive-bayes \
   --i-reference-reads $Outputs/database_seq.qza \
   --i-reference-taxonomy $Outputs/database_taxonomy.qza \
   --o-classifier $Outputs/taxo_classifier.qza \
   --verbose
 
# Assigning taxonomy with the trainned classifier
printf "\nAssigning taxonomy with the trainned classifier\n" 

qiime feature-classifier classify-sklearn \
  --i-classifier $Outputs/taxo_classifier.qza \
  --i-reads $Outputs/rep_seqs_dada2.qza \
  --o-classification $Outputs/taxonomy.qza

qiime metadata tabulate \
  --m-input-file $Outputs/taxonomy.qza \
  --o-visualization $Outputs/taxonomy.qzv

# Exporting tables and sequences
printf "\nExporting tables and sequences\n" 

qiime tools export \
  $Outputs/table_dada2.qza \
  --output-dir $Outputs/Data_analysis

biom convert -i $Outputs/Data_analysis/feature-table.biom  -o $Outputs/Data_analysis/feature-table.txt --to-tsv

qiime tools export \
   $Outputs/taxonomy.qza \
   --output-dir $Outputs/Data_analysis

qiime tools export \
  $Outputs/rep_seqs_dada2.qza \
  --output-dir $Outputs/rep_set2.fna
  
  
###################################################
###################################################
################ 18S data #########################
###################################################
###################################################

# Directories
Database_fasta='/home/olivierl/apps/Metagenetics_Projects/Databases/SILVA_132_QIIME_release/silva_132_99_18S.fna'
Database_txt='/home/olivierl/apps/Metagenetics_Projects/Databases/SILVA_132_QIIME_release/18S_taxonomy_7_levels.txt'
Metadata_file='/home/olivierl/lus/fish_farms/Metadata_18S.txt'
mkdir /home/olivierl/lus/fish_farms/Results_18S_Dada2
Outputs='/home/olivierl/lus/fish_farms/Results_18S_Dada2'
fastqManifest='/home/olivierl/lus/fish_farms/Fastq_manifest_18S.csv'

#Importing files into qiime format
printf "\nImporting fastq files to qiime2\n" 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $fastqManifest \
  --output-path $Outputs/paired_end_demux.qza \
  --source-format PairedEndFastqManifestPhred33

#Visualization of sequence quality
qiime demux summarize \
  --i-data $Outputs/paired_end_demux.qza \
  --o-visualization $Outputs/paired_end_demux.qzv \
  --verbose 
  
# Quality filtering, denoising and chimera filtering using DADA2
printf "\nStarting DADA2 process\n" 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $Outputs/paired_end_demux.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 225 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 216 \
  --p-chimera-method consensus \
  --p-min-fold-parent-over-abundance 1 \
  --p-n-threads 0 \
  --o-representative-sequences $Outputs/rep_seqs_dada2.qza \
  --o-table $Outputs/table_dada2.qza \
  --o-denoising-stats $Outputs/stats-dada2.qza \
  --verbose

printf "\nCreating visualisation files of DADA2 outputs\n" 
# Looking at the feature table and feature data summary
qiime feature-table summarize \
  --i-table $Outputs/table_dada2.qza \
  --o-visualization $Outputs/table_dada2.qzv \
  --m-sample-metadata-file $Metadata_file

qiime feature-table tabulate-seqs \
  --i-data $Outputs/rep_seqs_dada2.qza \
  --o-visualization $Outputs/rep_seqs_dada2.qzv
  
qiime metadata tabulate \
   --m-input-file $Outputs/stats-dada2.qza \
   --o-visualization $Outputs/stats-dada2.qzv
   
   
# Training a Naive Bayes classifier - assigning taxonomy from different databases
# Importing the fasta and taxonomy files
printf "\nInitiating taxonomic assignment\n" 
printf "\nImporting reference database\n" 

qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path $Database_fasta \
   --output-path $Outputs/database_seq.qza

 qiime tools import \
   --type 'FeatureData[Taxonomy]' \
   --source-format HeaderlessTSVTaxonomyFormat \
   --input-path $Database_txt \
   --output-path $Outputs/database_taxonomy.qza
   
# Extracting the sequence portion we have targeted
printf "\nExtracting the sequence portion we have targeted\n" 

 qiime feature-classifier extract-reads \
   --i-sequences $Outputs/database_seq.qza \
   --p-f-primer AGGGCAAKYCTGGTGCCAGC \
   --p-r-primer GRCGGTATCTRATCGYCTT \
   --o-reads $Outputs/database_seq.qza
   
# Trainning the classifier
printf "\nTrainning the classifier\n" 

 qiime feature-classifier fit-classifier-naive-bayes \
   --i-reference-reads $Outputs/database_seq.qza \
   --i-reference-taxonomy $Outputs/database_taxonomy.qza \
   --o-classifier $Outputs/taxo_classifier.qza \
   --verbose
 
# Assigning taxonomy with the trainned classifier
printf "\nAssigning taxonomy with the trainned classifier\n" 

qiime feature-classifier classify-sklearn \
  --i-classifier $Outputs/taxo_classifier.qza \
  --i-reads $Outputs/rep_seqs_dada2.qza \
  --o-classification $Outputs/taxonomy.qza

qiime metadata tabulate \
  --m-input-file $Outputs/taxonomy.qza \
  --o-visualization $Outputs/taxonomy.qzv

# Exporting tables and sequences
printf "\nExporting tables and sequences\n" 

qiime tools export \
  $Outputs/table_dada2.qza \
  --output-dir $Outputs/Data_analysis

biom convert -i $Outputs/Data_analysis/feature-table.biom  -o $Outputs/Data_analysis/feature-table.txt --to-tsv

qiime tools export \
   $Outputs/taxonomy.qza \
   --output-dir $Outputs/Data_analysis

qiime tools export \
  $Outputs/rep_seqs_dada2.qza \
  --output-dir $Outputs/rep_set2.fna
