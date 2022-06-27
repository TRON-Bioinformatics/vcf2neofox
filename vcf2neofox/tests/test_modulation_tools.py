from unittest import TestCase
import vcf2neofox.modulation_tools as modulation_tools


class TestModulationTools(TestCase):

    def test_build_neoantigen_no_mutation(self):
        neoantigen = modulation_tools.build_neoantigen('AAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
        self.assertIsNone(neoantigen)

    def test_build_neoantigen_mutation(self):
        neoantigen = modulation_tools.build_neoantigen(
            seq_wt_sequence='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            seq_mutated_sequence='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
        self.assertIsNotNone(neoantigen)
        self.assertIsNotNone(neoantigen.mutation.mutated_xmer)
        self.assertIsNotNone(neoantigen.mutation.wild_type_xmer)
        self.assertTrue(len(neoantigen.mutation.mutated_xmer) == len(neoantigen.mutation.wild_type_xmer))
        self.assertTrue(neoantigen.mutation.mutated_xmer != neoantigen.mutation.wild_type_xmer)

    def test_build_neoantigen_small_protein(self):
        neoantigen = modulation_tools.build_neoantigen(
            seq_wt_sequence='AAAAAAAAAAA',
            seq_mutated_sequence='AAAAADAAAAA')
        self.assertIsNotNone(neoantigen)
        self.assertEqual(neoantigen.mutation.mutated_xmer, 'AAAAADAAAAA')
        self.assertEqual(neoantigen.mutation.wild_type_xmer, 'AAAAAAAAAAA')
        self.assertTrue(len(neoantigen.mutation.mutated_xmer) == len(neoantigen.mutation.wild_type_xmer))
        self.assertTrue(neoantigen.mutation.mutated_xmer != neoantigen.mutation.wild_type_xmer)
