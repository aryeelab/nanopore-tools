module load local/singularity/2.5.2
export SINGULARITY_CACHEDIR=/data/aryee/sowmya/ctc_nanopore/images
(umask 002 && singularity pull docker://aryeelab/guppy)
rm -rf  ${SINGULARITY_CACHEDIR}/docker/*
rm -rf  ${SINGULARITY_CACHEDIR}/metadata/*
sleep 20s
(umask 002 && singularity pull docker://aryeelab/nanopore_minimap2)
rm -rf  ${SINGULARITY_CACHEDIR}/docker/*
rm -rf  ${SINGULARITY_CACHEDIR}/metadata/*
sleep 20s
(umask 002 && singularity pull docker://aryeelab/nanopolish:v0.11.0)
rm -rf  ${SINGULARITY_CACHEDIR}/docker/*
rm -rf  ${SINGULARITY_CACHEDIR}/metadata/*
sleep 20s
(umask 002 && singularity pull docker://aryeelab/nanopore_util)
rm -rf  ${SINGULARITY_CACHEDIR}/docker/*
rm -rf  ${SINGULARITY_CACHEDIR}/metadata/*
sleep 20s
