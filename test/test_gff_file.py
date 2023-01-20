import unittest
from module import gff_file


class MyTestCase(unittest.TestCase):

    def test_read_gff_file(self):
        mini_gff_file = './input_files/mini.gff'
        gene_organism = {'9519cfca-0c42-44d4-ab09-d37d33245d07': 'sand_box', '78fbaf52-0a8b-4acb-a93c-5d50238e47bf': 'sand_box'}

        mini_gff = gff_file.HandleGFF(mini_gff_file, gene_organism, '')
        mini_gff.read_gff_file()
        scaffold = mini_gff.fields[('scaffold', '9519cfca-0c42-44d4-ab09-d37d33245d07')]
        self.assertEqual(scaffold, '3R')

        strand = mini_gff.fields[('strand', '14f85617-72a1-4658-8abc-05f94e551114')]
        self.assertEqual(strand, '+')

        position = mini_gff.fields[('position', '26187b4d-cdb3-4263-a6af-a5556bfb5474_26')]
        self.assertEqual(position, (36675298, 36675597))

        gene_id = mini_gff.child_parent_relationship['5d7a93e5-6390-4568-832b-a5d7ca162ac6']
        self.assertEqual(gene_id, '9519cfca-0c42-44d4-ab09-d37d33245d07')

        mrna_parent_id = mini_gff.child_parent_relationship['14f85617-72a1-4658-8abc-05f94e551114']
        self.assertEqual('78fbaf52-0a8b-4acb-a93c-5d50238e47bf', mrna_parent_id)

        cds_parent_id = mini_gff.child_parent_relationship['79cdd16a-3988-4d9a-81a6-e0a17384013c']
        self.assertEqual('14f85617-72a1-4658-8abc-05f94e551114', cds_parent_id)

        gene_id_2 = mini_gff.get_gene_id('79cdd16a-3988-4d9a-81a6-e0a17384013c')
        self.assertEqual('78fbaf52-0a8b-4acb-a93c-5d50238e47bf', gene_id_2)

    def test_scan_true_gff_for_errors(self):
        true_gff_file = './input_files/true.gff'

        true_gene_organism = dict()

        with open('./input_files/two_gene_organism.tsv') as file_handle:
            for line in file_handle:
                gene_id, organism = line.rstrip().split("\t")
                true_gene_organism[gene_id] = organism

        true_gff = gff_file.HandleGFF(true_gff_file, true_gene_organism, '')
        true_gff.read_gff_file()

        true_gff.scan_gff_for_errors()
        self.assertEqual(true_gff.errors, {})

    def test_scan_false_gff_for_errors(self):
        false_gff_file = './input_files/simple_false.gff'

        false_gene_organism = dict()
        with open('./input_files/simple_organism.tsv') as file_handle:
            for line in file_handle:
                gene_id, organism = line.rstrip().split("\t")
                false_gene_organism[gene_id] = organism

        false_gff = gff_file.HandleGFF(false_gff_file, false_gene_organism, '')
        false_gff.read_gff_file()

        false_gff.scan_gff_for_errors()

        self.assertEqual(4, len(false_gff.errors))
        scaffold_error = false_gff.errors['f18e9140-5589-405b-86c3-0e54cf01390d']
        self.assertEqual(1, len(scaffold_error.gff_format_error))
        self.assertEqual('simple@ebi.ac.uk', scaffold_error.owner)
        self.assertEqual('sand_box', scaffold_error.organism_name)
        self.assertEqual('6c797cb7-1246-4add-aab4-a3287a12d27a', scaffold_error.gene_id)
        self.assertEqual('8c5922b3-fe26-49dc-8aeb-9ac867b96c07', scaffold_error.mrna_id)
        self.assertEqual('AGAP010269', scaffold_error.gene_name)
        self.assertEqual('3R:51887721..51894443', scaffold_error.locus)
        self.assertEqual('scaffold', scaffold_error.gff_format_error[0].get('field_type'))
        self.assertEqual('exon', scaffold_error.gff_format_error[0].get('feature_type'))
        self.assertEqual('f18e9140-5589-405b-86c3-0e54cf01390d', scaffold_error.gff_format_error[0].get('feature_id'))
        self.assertEqual('2R', scaffold_error.gff_format_error[0].get('feature_value'))
        self.assertEqual('8c5922b3-fe26-49dc-8aeb-9ac867b96c07', scaffold_error.gff_format_error[0].get('parent_id'))
        self.assertEqual('3R', scaffold_error.gff_format_error[0].get('parent_value'))

if __name__ == '__main__':
    unittest.main()
