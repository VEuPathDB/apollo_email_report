import unittest
from module import transcript


class MyTestCase(unittest.TestCase):
    # There needs to be a new test for webservice endpoints.
    # def test_has_sequence(self):
        # mrna = transcript.CodingSequence('f3687329-5172-48d1-bce7-134d914f1d28', 'sandbox_arabiensis', 'SDBKB704125')
        # mrna.get_sequence('https://apollo.vectorbase.org/Apollo/')
        # self.assertEqual('ATGTCA', mrna.sequence[0:6])

    def test_has_start_codon(self):
        mrna = transcript.CodingSequence('b0b85443-0ff0-4fec-b919-f7bbeb626072', 'sandbox_arabiensis', 'SDBKB704125')
        mrna.sequence = "ATGCACTGA"
        mrna.coding_sequence_has_start_codon()
        self.assertEqual(True, mrna.coding_sequence_has_start_codon())

        mrna.sequence = "ATCCACTGA"
        mrna.coding_sequence_has_start_codon()
        self.assertEqual('no start codon', mrna.errors['start_codon'])

    def test_has_stop_codon(self):
        mrna = transcript.CodingSequence('b0b85443-0ff0-4fec-b919-f7bbeb626072', 'sandbox_arabiensis', 'SDBKB704125')

        mrna.sequence = "ATGCACTGA"  # TGA amber stop codon
        self.assertEqual(True, mrna.coding_sequence_has_stop_codon())

        mrna.sequence = "ATGCAGTAA"  # TAA ochre stop codon
        self.assertEqual(mrna.coding_sequence_has_stop_codon(), True)

        mrna.sequence = "ATGCAGTAG"  # TAG opal stop codon
        self.assertEqual(mrna.coding_sequence_has_stop_codon(), True)

        mrna.sequence = "ATGCAGCAG"  # NO stop codon
        self.assertEqual(mrna.coding_sequence_has_stop_codon(), False)
        self.assertEqual('no stop codon', mrna.errors['stop_codon'])

    def test_no_internal_stop_codon(self):
        mrna = transcript.CodingSequence('b0b85443-0ff0-4fec-b919-f7bbeb626072', 'sandbox_arabiensis', 'SDBKB704125')

        mrna.sequence = "ATGCACCTCGAGTAA"
        mrna.coding_sequence_no_internal_stop_codon()
        self.assertEqual(mrna.coding_sequence_no_internal_stop_codon(), True)

        mrna.sequence = "ATGCACTAACTCGAGTAA"  # TAA
        mrna.coding_sequence_no_internal_stop_codon()
        self.assertEqual(mrna.errors['no_internal_stop_codon'], '1 internal stop codon')

        mrna.sequence = "ATGCACCTCGAGTAGTAA"  # TAG
        mrna.coding_sequence_no_internal_stop_codon()
        self.assertEqual(mrna.errors['no_internal_stop_codon'], '1 internal stop codon')

        mrna.sequence = "ATGTGACACCTCTAATAA"  # TGA,TAA
        mrna.coding_sequence_no_internal_stop_codon()
        self.assertEqual(mrna.errors['no_internal_stop_codon'], '2 internal stop codon')


if __name__ == '__main__':
    unittest.main()
