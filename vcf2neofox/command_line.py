import argparse

import neofox.model.conversion
import pandas as pd
import vcf2neofox
from logzero import logger
import vcf2neofox.modulation_tools as tools
import vcf2neofox.data_loading as data_loading
from cyvcf2 import VCF

epilog = "Copyright (c) 2022 TRON gGmbH (See LICENSE for licensing details)"


def vcf2neofox_cli():
    # set up logger
    parser = argparse.ArgumentParser(description="vcf2neofox v{}".format(vcf2neofox.VERSION),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=epilog)
    parser.add_argument("--output-table", dest="output_table", action="store",
                        help="The path to the output table in NeoFox format", required=True)
    parser.add_argument(
        "--patient-id",
        dest="patient_id",
        help="Patient identifier",
        required=True,
    )
    parser.add_argument(
        "--vcf",
        dest="vcf",
        help="Mutations",
        required=True,
    )
    parser.add_argument(
        "--fasta",
        dest="fasta",
        help="reference genome",
        required=True,
    )
    parser.add_argument(
        "--gtf",
        dest="gtf",
        help="gene annotations",
        required=True,
    )
    args = parser.parse_args()
    
    logger.info("VCF2neofox starting...")
    
    # VCF Input File
    vcf_file = args.vcf
    
    # Read references
    logger.info("Reading references...")
    reference_genome, gtf = data_loading.load_references(gtf=args.gtf, fasta=args.fasta)
    
    # Load VCF
    logger.info("Reading VCF...")
    vcf = VCF(vcf_file)
    neoantigens = []
    # Iterate variants
    for v in vcf:
        
        variant = v
        # get overlapping transcripts
        transcripts = tools.get_overlapping_transcripts(gtf, variant)
        # Stop if the variant is out of transcript regions
        if transcripts.empty:
            continue
        
        # Check overlapping exons for each transcript
        for tr_id in transcripts.transcript_id:
            exons = tools.get_exons(gtf, reference_genome, tr_id)
            if (not tools.is_coding(exons, variant)) | (exons is None):
                continue
            
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
            gene = exons.gene_name.drop_duplicates().iloc[0]
            
            # Create neoantigen model
            try:
                neoantigen = tools.build_neoantigen(normal_protein, mutated_protein, patient_identifier=args.patient_id, gene = gene)
            except Exception as ex:
                # skip failing Neoantigen (eg: sequences containing * or invalid AAs)
                logger.warning(ex)
                continue
            
            # Create output table
            if (neoantigen is not None) & (neoantigen not in neoantigens):
                neoantigens.append(neoantigen)

    # writes a CSV with neoantigens
    neoantigens_df = pd.json_normalize(
        data=[n.to_dict(include_default_values=True) for n in neoantigens]
    )
    neoantigens_df.to_csv(args.output_table, sep='\t', index=False)

    logger.info("VCF2neofox finished!")
