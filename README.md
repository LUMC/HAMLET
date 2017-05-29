# Hamlet

Hamlet is a pipeline for analysis of human acute myeloid leukemia RNA-seq samples.


## Running Hamlet

    snakemake -p -T --use-conda -s hamlet.snake --configfile config.yml \
        --rerun-incomplete --restart-times 3 \
        --cluster-config config-cluster.yml --jobname 'hamlet.{jobid}' \
        --jobs 100 -w 120 --max-jobs-per-second 3 \
        --drmaa ' -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -cwd -V' \
        --drmaa-log-dir .drmaa-logs
