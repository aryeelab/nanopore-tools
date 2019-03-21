run_id=test-run-1
fast5_zip=test_data/test-run-1.zip
flowcell=FLO-MIN106      
kit=SQK-RBK004


# Unzip fast5
fast5_path=`jar tf $fast5_zip | grep 'fast5/$'` # find path to fast5 dir within zip
echo "Unzipping $fast5_zip"
jar -xf $fast5_zip
echo "Moving $fast5_path to `basename $fast5_path`"
mv $fast5_path .


# Guppy basecall and demultiplex
guppy_basecaller -r -i /fast5 -s guppy_basecaller -q 0 --flowcell $flowcell --kit $kit --cpu_threads_per_caller 2 --qscore_filtering --min_qscore 7
guppy_barcoder -i guppy_basecaller/pass -s guppy_barcoder --barcode_kits ${kit} 