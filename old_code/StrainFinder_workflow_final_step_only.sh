#!/bin/bash

# Only argument should be path to config file that holds all the relevant settings.

. $1

module load python/2
source $HOME/.strainfinder-virtualenv/bin/activate

cd $PARENT_WORKING_DIR/$WORKING_DIR_NAME

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
                                                                 --em em.{}.cpickle \
                                                                 --em_out em.{}.cpickle \
                                                                 --otu_out otu_table.{}.txt \
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

