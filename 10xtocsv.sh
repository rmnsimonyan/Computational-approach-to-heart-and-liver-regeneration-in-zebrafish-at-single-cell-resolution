# Copy the script into the directory of data and run from there
dir=$1
cp 10xGenmtx_and_tsv_to_csv.sh $dir
cd $dir
sh 10xGenmtx_and_tsv_to_csv.sh
