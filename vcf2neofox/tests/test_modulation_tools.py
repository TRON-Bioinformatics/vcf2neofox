from unittest import TestCase
import vcf2neofox.modulation_tools as modulation_tools


class TestModulationTools(TestCase):

    def test_build_neoantigen_no_mutation(self):
        neoantigen = modulation_tools.build_neoantigen(
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA', patient_identifier='123'
        )
        self.assertIsNone(neoantigen)

    def test_build_neoantigen_mutation(self):
        neoantigen = modulation_tools.build_neoantigen(
            seq_wt_sequence='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            seq_mutated_sequence='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            patient_identifier='123')
        self.assertIsNotNone(neoantigen)
        self.assertIsNotNone(neoantigen.mutated_xmer)
        self.assertIsNotNone(neoantigen.wild_type_xmer)
        self.assertTrue(len(neoantigen.mutated_xmer) == len(neoantigen.wild_type_xmer))
        self.assertTrue(neoantigen.mutated_xmer != neoantigen.wild_type_xmer)

    def test_build_neoantigen_small_protein(self):
        neoantigen = modulation_tools.build_neoantigen(
            seq_wt_sequence='AAAAAAAAAAA',
            seq_mutated_sequence='AAAAADAAAAA',
            patient_identifier='123')
        self.assertIsNotNone(neoantigen)
        self.assertEqual(neoantigen.mutated_xmer, 'AAAAADAAAAA')
        self.assertEqual(neoantigen.wild_type_xmer, 'AAAAAAAAAAA')
        self.assertTrue(len(neoantigen.mutated_xmer) == len(neoantigen.wild_type_xmer))
        self.assertTrue(neoantigen.mutated_xmer != neoantigen.wild_type_xmer)
