# Create a subset of the arriba database, which overlaps the hamlet-ref.bed file
bed=test/data/reference/hamlet-ref.bed

# We don't care about the blacklist for now, so we just take a few lines
zcat blacklist_hg38_GRCh38_v2.4.0.tsv.gz | head | gzip > test/blacklist_hg38_GRCh38_v2.4.0.head.tsv.gz

# Get known fusions, BCR only
zgrep -A 1 BCR known_fusions_hg38_GRCh38_v2.4.0.tsv.gz | gzip > test/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz

# Get the relevant protein domains that overlap hamlet-ref.bed
# First, append 'chr' to all the chromosome names
awk '{ print "chr" $0}' protein_domains_hg38_GRCh38_v2.4.0.gff3 > protein_domains_hg38_GRCh38_v2.4.0.chr.gff3

# Next, use bedtools to intersect
bedtools intersect -a protein_domains_hg38_GRCh38_v2.4.0.chr.gff3 -b $bed > protein_domains_hg38_GRCh38_v2.4.0.chr.isect.gff3

# Cut off the 'chr' prefix again
awk '{ print substr($0,4,length($0))}' protein_domains_hg38_GRCh38_v2.4.0.chr.isect.gff3 > test/protein_domains_hg38_GRCh38_v2.4.0.gff3
