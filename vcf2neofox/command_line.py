import argparse
import vcf2neofox
import logging

epilog = "Copyright (c) 2022 TRON gGmbH (See LICENSE for licensing details)"


def vcf2neofox_cli():
    # set up logger
    parser = argparse.ArgumentParser(description="vcf2neofox v{}".format(vcf2neofox.VERSION),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=epilog)
    # TODO: add arguments
    args = parser.parse_args()

    logging.info("VCF2neofox starting...")
    # TODO: call your code
    logging.info("VCF2neofox finished!")