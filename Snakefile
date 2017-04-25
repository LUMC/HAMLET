from rattle import Run

RUN = Run(config)

include: "includes/qc/Snakefile"

rule all:
    input:
        fqs=[RUN.output("{sample}/{sample}-{pair}.fq.gz", fmt=True,
                        sample=unit.sample, pair=pair)
             for unit in RUN.unit_names for pair in ("R1", "R2")],
