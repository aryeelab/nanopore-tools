
    gcloud compute instances create megalodon \
        --machine-type a2-highgpu-2g \
        --zone us-central1-a \
        --boot-disk-size 2000GB \
        --image-family cos-97-lts \
        --image-project cos-cloud \
        --maintenance-policy TERMINATE --restart-on-failure 
        
        gcloud compute ssh --zone us-central1-a megalodon
        
        # Localize inputs
        # Reads
        toolbox gsutil cp gs://aryeelab-nanopore/large-test/HEPG2-H3K27me3test.tar.gz /media/root/home/$USER/
        # Genome FASTA
        toolbox gsutil cp gs://aryeelab-nanopore/griffinlab/references/GCA_000001405.15_GRCh38_full_analysis_set.fna /media/root/home/$USER/
        # Models
        toolbox gsutil cp gs://aryeelab-nanopore/griffinlab/references/basecall.tar.gz /media/root/home/$USER/
        
        # Install GPU drivers
        sudo cos-extensions install gpu
        sudo mount --bind /var/lib/nvidia /var/lib/nvidia
        sudo mount -o remount,exec /var/lib/nvidia

        # Configure artifact registry credentials
        docker-credential-gcr configure-docker --registries us-central1-docker.pkg.dev

        docker run --rm -it \
          --volume /home/$USER:/work \
          --volume /var/lib/nvidia/lib64:/usr/local/nvidia/lib64 \
          --volume /var/lib/nvidia/bin:/usr/local/nvidia/bin \
          --device /dev/nvidia0:/dev/nvidia0 \
          --device /dev/nvidia-uvm:/dev/nvidia-uvm \
          --device /dev/nvidiactl:/dev/nvidiactl \
          us-central1-docker.pkg.dev/aryeelab/docker/megalodon
          
		# In container
        cd /work  
          
        reads="HEPG2-H3K27me3test.tar.gz"
        model="basecall.tar.gz"
        genome="GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"
        modmotifs="--mod-motif Z CG 0 --mod-motif Y A 0"
        device="cuda:all --processes 10"

        nvidia-smi
        mkdir ./in
        tar zxvf ${reads} -C ./in
        mkdir ./out
        mkdir ./basecall_models
        tar zxvf ${model} -C ./basecall_models
        time megalodon \
        ./in \
        --output-directory "./out" \
        --overwrite \
        --guppy-server-path /usr/bin/guppy_basecall_server \
        --guppy-params "-d ./basecall_models" \
        --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
        --outputs basecalls mappings mod_mappings mods per_read_mods \
        --reference ${genome} ${modmotifs} \
        --devices ${device}
        
        