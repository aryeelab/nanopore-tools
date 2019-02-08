sample_id=test-sample-1
fast5_zip=test_data/test-sample-1.zip

flowcell=FLO-MIN106      
kit=SQK-RBK004


# Unzip fast5
fast5_path=`jar tf $fast5_zip | grep 'fast5/$'` # find path to fast5 dir within zip
jar -xf $fast5_zip
mv $fast5_path .


# Albacore
source /apps/lab/aryee/pyenv/versions/nanopore/bin/activate
read_fast5_basecaller.py -f  $flowcell -k $kit -t 1 -r -i  fast5 -o fastq -q 0 -s albacore


  

