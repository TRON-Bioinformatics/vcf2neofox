from pysam import FastaFile
import gtfparse

# Reference assemblies
# HG38: /projects/data/human/ensembl/GRCh38.97/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# HG37: /projects/data/human/ensembl/GRCh37.90/Homo_sapiens.GRCh37.dna.primary_assembly.fa
fasta_file_grch37 = "data/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
fasta_file_grch38 = "data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# GTF files
# HG38: /projects/data/human/ensembl/GRCh38.97/Homo_sapiens.GRCh38.97.gtf
# HG37: /projects/data/human/ensembl/GRCh37.90/Homo_sapiens.GRCh37.87.gtf --> released with dna GRCh37.90
gtf_file_grch37 = "data/Homo_sapiens.GRCh37.87.gtf"
gtf_file_grch38 = "data/Homo_sapiens.GRCh38.97.gtf"

# TODO: Only read selected columns to make things faster


def load_references(fasta, gtf):
    reference_genome = FastaFile(fasta)
    gtf = gtfparse.read_gtf(gtf)
    filtered_gtf = gtf[(gtf.gene_biotype == "protein_coding") & (gtf.tag.str.contains("basic"))]
    return reference_genome, filtered_gtf
