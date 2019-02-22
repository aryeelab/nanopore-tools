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


# Albacore basecall and demultiplex
source /apps/lab/aryee/pyenv/versions/nanopore/bin/activate
read_fast5_basecaller.py -f  $flowcell -k $kit -t 1 -r -i  fast5 -o fastq -q 0 -s albacore


  

