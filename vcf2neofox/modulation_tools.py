from Bio.Seq import Seq
from cyvcf2 import VCF
from neofox.model.factories import NeoantigenFactory
from neofox.model.neoantigen import Neoantigen


def get_overlapping_transcripts(gtf, VCF_item):
    # HG19: GTF does not contain "chr" in the chromosome field, HG38 does
    # TODO: Is the trancsript status relevant --> "putative"; highest transcript level?
    overlapping_transcripts = gtf[(gtf.seqname == VCF_item.CHROM.strip("chr")) & 
                                (gtf.start <= VCF_item.POS) &
                                (gtf.end > VCF_item.POS) & 
                                (gtf.feature == "transcript")]
    return(overlapping_transcripts[['seqname',  'start', 'end', "gene_name", "transcript_id", "strand"]])


# fetch all exons for a given transcript
# NOTE: beware to read elements of feature_type=CDS and not feature_type=exon, the later include the UTR regions which we do not want
def get_exons(gtf, reference_genome, transcript_id):
    exons = gtf[(gtf.transcript_id == transcript_id) & (gtf.feature == "CDS")].sort_values("start")
    exon_sequences = []
    for _, exon in exons.iterrows():
        # NOTE: beware that coordinates in the GTF are 1-based but when querying the genome when pysam we need 0-based
        exon_sequences.append(reference_genome.fetch(exon.seqname.strip("chr"), exon.start - 1, exon.end))    
    exons["sequence"] = exon_sequences
    exons["sequence_len"] = exons["sequence"].transform(len)
    exons = exons[["gene_name", "exon_id", "seqname", "start", "end", "strand", "exon_number", "sequence", "sequence_len"]]
    return(exons)


# Check if the variant is in the coding region
def is_coding(exons, variant):
    overlapping_exon = exons[(exons.seqname == variant.CHROM.strip("chr")) & (exons.start <= variant.POS) & (exons.end > variant.POS)]
    if overlapping_exon.empty:
        return False
    else:
        return True


def create_coding_regions(exons, variant):
    # Create original DNA
    dna_sequence = ""
    for _, exon in exons.iterrows():
        dna_sequence = dna_sequence + exon.sequence
    
    # Altered transcript
    overlapping_exon = exons[(exons.seqname == variant.CHROM.strip("chr")) & (exons.start <= variant.POS) & (exons.end > variant.POS)]
    overlapping_sequence = overlapping_exon.sequence.iloc[0]
    variant_relative_position = variant.POS - overlapping_exon.start.iloc[0]
    mutated_sequence=overlapping_sequence[0: variant_relative_position] + variant.ALT[0].lower() + overlapping_sequence[variant_relative_position + 1 : ]
    
    # Create mutated DNA
    mutated_dna_sequence = ""
    for _, exon in exons.iterrows():
        if exon.exon_number == overlapping_exon.exon_number.iloc[0]:
            mutated_dna_sequence = mutated_dna_sequence + mutated_sequence
        else:
            mutated_dna_sequence = mutated_dna_sequence + exon.sequence
    
    return(dna_sequence, mutated_dna_sequence)


# Translate DNA sequences to protein
def translate(dna, strand):
    protein = None
    if strand == "+":
        protein = Seq(dna).translate()
    elif strand == "-":
        protein = Seq(dna).reverse_complement().translate()
    else:
        print("Wrong strand argument!")
    return(protein)


# Generate epitope sequence
def build_neoantigen(seq_mutated_sequence: str, seq_wt_sequence: str, patient_identifier: str = "p", gene: str = "") -> Neoantigen:
    # Find the mutated aminoacid in the protein and generate the neoepitopes
    aa_position = 0
    for ref, alt in zip(seq_wt_sequence, seq_mutated_sequence):
        if ref != alt:
            break
        aa_position += 1
    
    # Skip synonymous variants
    neoantigen = None
    if(aa_position != len(seq_wt_sequence)):
        mutated_xmer = seq_mutated_sequence[max(aa_position - 13, 0): min(aa_position + 14, len(seq_wt_sequence))]
        wt_xmer = seq_wt_sequence[max(aa_position - 13, 0): min(aa_position + 14, len(seq_wt_sequence))]
        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer=mutated_xmer,
            wild_type_xmer=wt_xmer,
            patient_identifier=patient_identifier,
            gene=gene,
            position=14
        )

    return neoantigen
