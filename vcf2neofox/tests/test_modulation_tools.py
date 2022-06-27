from unittest import TestCase
import vcf2neofox.modulation_tools as modulation_tools


class TestModulationTools(TestCase):

    def test_no_mutation(self):
        neoantigen = modulation_tools.build_neoantigen('AAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
        self.assertIsNone(neoantigen)

