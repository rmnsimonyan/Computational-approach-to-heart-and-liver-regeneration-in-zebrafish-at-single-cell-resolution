# Copy the script into the directory of data and run from there
dir=$1
cd data
unzip $dir.zip
cd ..
cp 10xGenmtx_and_tsv_to_csv.sh data/$dir
cd data/$dir
sh 10xGenmtx_and_tsv_to_csv.sh
