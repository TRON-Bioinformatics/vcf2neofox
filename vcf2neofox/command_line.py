import argparse

import neofox.model.conversion

import vcf2neofox
import logging
import vcf2neofox.modulation_tools as tools
import vcf2neofox.data_loading as data_loading
from cyvcf2 import VCF
from neofox.model.neoantigen import Neoantigen, Mutation

epilog = "Copyright (c) 2022 TRON gGmbH (See LICENSE for licensing details)"

def vcf2neofox_cli():
    # set up logger
    parser = argparse.ArgumentParser(description="vcf2neofox v{}".format(vcf2neofox.VERSION),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=epilog)
    # TODO: add arguments
    args = parser.parse_args()
    
    logging.info("VCF2neofox starting...")
    
    # VCF Input File
    vcf_file = "data/H018-RM7N.varscan2.hc.vcf"
    
    # Read references
    reference_genome, gtf = data_loading.load_references('grch37')
    
    # Load VCF
    vcf = VCF(vcf_file)
    neoantigens = []
    # Iterate variants
    for v in vcf:
        
        variant = v
        transcripts = None
        # get overlapping transcripts
        transcripts = tools.get_overlapping_transcripts(gtf, variant)
        # Stop if the variant is out of transcript regions
        if transcripts.empty:
            continue
        
        # Check overlapping exons for each transcript
        for tr_id in transcripts.transcript_id:
            exons = None
            exons = tools.get_exons(gtf, reference_genome, tr_id)
            if (not tools.is_coding(exons, variant)) | (exons is None):
                continue
            
            # Create dna sequence for normal and mutated protein
            dna_sequence, mutated_sequence = tools.create_coding_regions(exons, variant)
            strand = None
            if exons.strand.unique() == "-":
                strand = "-"
            if exons.strand.unique() == "+":
                strand = "+"
            
            normal_protein = tools.translate(dna_sequence, strand)
            mutated_protein = tools.translate(mutated_sequence, strand)
            
            neoantigen = tools.build_neoantigen(normal_protein, mutated_protein)
            
            # Create output table
            if neoantigen is not None:
                neoantigens.append(neoantigen)

    # writes a CSV with neoantigens
    neoantigens_df = neofox.model.conversion.ModelConverter().annotations2table(neoantigens)
    neoantigens_df.to_csv('your_neoantigens', sep='\t', index=False)

    logging.info("VCF2neofox finished!")