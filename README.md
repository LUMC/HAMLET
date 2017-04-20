# Hamlet

Hamlet is a pipeline for analysis of human acute myeloid leukemia RNA-seq samples.


## Running Hamlet

    snakemake --use-conda -p -s qc.snake --configfile config.yml-j 2 --cores 4 --jobname 'hamlet.{jobid}' --drmaa ' -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -cwd -V' -T --rerun-incomplete --jobs 100 -w 120 --max-jobs-per-second 3 --cluster-config config-cluster.yml
