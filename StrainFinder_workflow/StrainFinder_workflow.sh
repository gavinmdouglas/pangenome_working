#!/bin/bash

# Only argument should be path to config file that holds all the relevant settings.

. $1

mkdir -p $PARENT_WORKING_DIR/$WORKING_DIR_NAME

cd $PARENT_WORKING_DIR/$WORKING_DIR_NAME

module load python/3

# Specify which samples to retain for this analysis
# (e.g., only those where the genes were called as present by pandora)
python $PATH_TO_SCRIPTS/identify_samples_w_genes_present.py \
	-p $PRESENCE_ABSENCE_MATRIX \
	-g $INPUT_GENE_IDS \
	-o samples_to_keep.txt

# Get FASTA of genes of interest
python $PATH_TO_SCRIPTS/gene_fasta_to_gene_format_for_kpileup.py \
       -f $MASTER_FASTA \
       -g $INPUT_GENE_IDS \
       -o gene_kpileup_infiles

# Do this step in case multiple genes specified.
cat gene_kpileup_infiles/*gene > input.gene
rm -r gene_kpileup_infiles

mkdir kpileup_out
for SAMPLE in $( cat samples_to_keep.txt ); do
    BAMFILE="$BAM_FILE_DIR/$SAMPLE.$BAM_SUFFIX"
    echo "perl $STRAINFINDER_FOLDER/preprocess/3.kpileup.pl $SAMPLE $BAMFILE input.gene 20 0 $MIN_DEPTH > kpileup_out/$SAMPLE.kpileup" >> kpileup_cmds.sh
done

cat kpileup_cmds.sh | parallel -j $NUM_CORE '{}'

### source /home/gdouglas/local/miniconda3/etc/profile.d/conda.sh
### conda activate strainfinder

module load python/2
source $HOME/.strainfinder-virtualenv/bin/activate

python $STRAINFINDER_FOLDER/preprocess/4.kp2np_edit.py \
       --samples samples_to_keep.txt \
       --gene_file input.gene \
       --out alignments.cPickle

for ID in $( cat $INPUT_GENE_IDS ); do echo -e "$WORKING_DIR_NAME\t$ID"; done > mock_genome2contig.tsv

python $STRAINFINDER_FOLDER/preprocess/5.filter_np_EDIT.py --aln alignments.cPickle \
                                                           --map mock_genome2contig.tsv \
                                                           --samples samples_to_keep.txt \
                                                           --tlen 0


mkdir strain_fitting
cd strain_fitting

FILTERED_ALIGN="$WORKING_DIR_NAME.np.cPickle"

for i in $( seq 2 $MAX_NUM_STRAINS_TO_TEST); do
	echo "python $STRAINFINDER_FOLDER/StrainFinder_edit.py \
                                                                 --aln ../$FILTERED_ALIGN \
                                                                 -N $i \
                                                                 --max_reps $MAX_REPS \
                                                                 --dtol $DTOL \
                                                                 --ntol $NTOL \
                                                                 --max_time $MAX_TIME \
                                                                 --converge \
                                                                 --em em.$i.cpickle \
                                                                 --em_out em.$i.cpickle \
                                                                 --otu_out otu_table.$i.txt \
                                                                 --log log.txt \
                                                                 --n_keep $N_KEEP \
                                                                 --force_update \
                                                                 --merge_out \
                                                                 --msg" >> ../StrainFinder_cmds.sh
done


cat ../StrainFinder_cmds.sh | parallel -j $NUM_CORE '{}'

python $PATH_TO_SCRIPTS/summarize_strainfinder_AICs.py \
             -i strain_fitting \
             -o strain_fit_summary.tsv

