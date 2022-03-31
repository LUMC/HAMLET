# Reference data source for STAR-fusion
CTAT := GRCh38_gencode_v22_CTAT_lib_Mar012021.source
CTAT_URL := https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/$(CTAT).tar.gz
REF_chrM := ../../test/data/reference/chrM.fa

# Gencode file, which is part of the STAR-fusion reference data
GENCODE_VERSION := v22
GENCODE := $(CTAT)/gencode.$(GENCODE_VERSION).annotation.gtf
GENCODE_chrM := $(CTAT)/gencode.$(GENCODE_VERSION).annotation.chrM.gtf

# fusion_lib file, which contains known fusion events.
FUSION_LIB := $(CTAT)/fusion_lib.Mar2021.dat.gz
FUSION_LIB_chrM := $(CTAT)/fusion_lib.Mar2021.dat.chrM.gz

# STAR-fusion singularity image to use (the version from biocontainers is
# missing some tools that are needed to build the reference set)
STAR := star-fusion.v1.10.1.simg
STAR_URL := https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/$(STAR)

# Tools
# (located inside the singularity container)
PREP_GENOME := /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl

# Output files
GENOME_LIB := ctat_genome_lib_build_dir_chrM
STAR_INDEX := $(GENOME_LIB)/ref_genome.fa.star.idx

#
# Maintenance targets.
#

.PHONY: all clean

all: $(CTAT) $(STAR) $(GENCODE_chrM) $(FUSION_LIB_chrM) $(GENOME_LIB) $(STAR_INDEX)

clean:
	rm -rf ctat_genome_lib_build_dir_chrM/ _* *.hmm.* *.hmm ref_annot.* $(CTAT) $(CTAT).tar.gz PFAM.domtblout.dat.gz $(STAR) pipeliner.*


#
# Reference data download targets.
#

$(CTAT):
	wget -c $(CTAT_URL) && tar -xvf $@.tar.gz

$(STAR):
	wget -c $(STAR_URL)

# 
# Subset there reference files to limit them to chrM
#

# Subset the gencode file to only contain chrM
$(GENCODE_chrM): $(GENCODE)
	grep "^chrM" $< > $@

# Subset the fusion lib file to only contain fusions between two genes that are
# both located on chrM
$(FUSION_LIB_chrM): $(GENCODE_chrM) $(FUSION_LIB)
	python fusion-lib-subset.py --gtf $(word 1, $^) --input $(word 2, $^) --output $@

#
# Build the required reference files
#

# Build the genome library for STAR-fusion
$(GENOME_LIB): $(STAR) $(REF_chrM) $(GENCODE_chrM) $(FUSION_LIB_chrM)
	singularity exec $(word 1, $^) \
		$(PREP_GENOME) \
		--fusion_annot_lib CTAT \
		--pfam_db current \
		--dfam_db human \
		--genome_fa $(word 2, $^) \
		--gtf $(word 3, $^) \
		--fusion_annot_lib $(word 4, $^) \
		--output_dir $@
	
	# Remove the STAR index, since we need to re-run that with modified settings
	# for small genomes
	rm -rf $@/ref_genome.fa.star.idx

# Build the STAR index, for small genomes
$(STAR_INDEX): $(STAR) $(REF_chrM) $(GENCODE_chrM)
	singularity run $(word 1, $^) \
		STAR \
		--runMode genomeGenerate \
		--limitGenomeGenerateRAM 40419136213 \
		--genomeChrBinNbits 16 \
		--genomeSAindexNbases 6 \
	    --sjdbOverhang 150 \
		--genomeDir $@ \
		--genomeFastaFiles $(word 2, $^) \
		--sjdbGTFfile $(word 3, $^)
