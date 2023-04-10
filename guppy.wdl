workflow ont_guppy_gpu {

  	File fast5_tar
  	String basecall_ref
  	Int disksize = 375
  
	call basecall {
  	input:
    	fast5_tar = fast5_tar,
        basecall_ref = basecall_ref,
        disksize=disksize
    }
	output {
    	File output_bases = basecall.output_bases
    }
}    
        
task basecall { 
	String basecall_ref
    File fast5_tar
    Int disksize
    
    command {
        set -o errexit
        set -o nounset
        set -o pipefail

        tar -xzvf ${fast5_tar}
        guppy_basecaller -i fast5 -s output_bases -c ${basecall_ref} -x "cuda:0" --num_callers 14 --gpu_runners_per_device 8 --compress_fastq
        tar -czvf output_bases.tar.gz output_bases
  	}
  
  	runtime {
      continueOnReturnCode: false
      docker: "quay.io/kdong2395/guppy-gpu:v0.2"
      disks: "local-disk ${disksize} LOCAL"
      bootDiskSizeGb: 25
      gpuType: "nvidia-tesla-p100"
      gpuCount: 1
      memory: "64GB"
      cpu: 16
      zones: "us-central1-c"
  	}
	output {
      File output_bases = "output_bases.tar.gz"
  	}
}
