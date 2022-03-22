
GENE=$1

GENE_SPLIT=($(echo $GENE | tr "_" "\n"))

SPECIES=${GENE_SPLIT[0]}"_"${GENE_SPLIT[1]}

ALN_FILE="/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/prepped_input/"$GENE".np.cPickle"

mkdir -p /data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_output/$SPECIES/$GENE/strain_running

cd /data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_output/$SPECIES/$GENE/strain_running

for NUM in {1..30}; do

    python /home/gdouglas/local/prg/StrainFinder/StrainFinder_edit.py \
                                                             --aln $ALN_FILE \
                                                             -N $NUM \
                                                             --max_reps 100 \
                                                             --dtol 1 \
                                                             --ntol 2 \
                                                             --max_time 60 \
                                                             --converge \
                                                             --em em.$GENE.$NUM.cpickle \
                                                             --em_out em.$GENE.$NUM.cpickle \
                                                             --otu_out otu_table.$GENE.$NUM.txt \
                                                             --log log.$GENE.$NUM.txt \
                                                             --n_keep 3 \
                                                             --force_update \
                                                             --merge_out \
                                                             --msg
done

cd ..

SUMMARY_OUT=$GENE".strain_fit_summary.tsv"
python /home/gdouglas/scripts/pangenome_working/StrainFinder_workflow/summarize_strainfinder_AICs.py -i $PWD/strain_running/ \
                                                                                                     -o $SUMMARY_OUT

CHECK_OUT=$GENE".check.tsv"
python /home/gdouglas/scripts/pangenome_working/StrainFinder_workflow/check_StrainFinder_output_dir.py -g $GENE \
                                                                                                       --parentfolder /data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_output/$SPECIES/ \
                                                                                                       --subfolder strain_running/ \
                                                                                                       --max_strains 30 > $CHECK_OUT

OTU_FIT=$( grep "passed" $GENE.check.tsv | awk '{ print $4 }' )

if [ ! -z $OTU_FIT ]; then
    cp $OTU_FIT ./;
fi

tar -zcvf strain_running.tar.gz strain_running --remove-files

