import unittest
import filecmp
import datetime
from module import annotation_quality_report, gff_file


class MyTestCase(unittest.TestCase):

    def test_writing_out_error(self):
        false_gff_file = './input_files/simple_false.gff'

        false_gene_organism = dict()
        file_handle = open('./input_files/simple_organism.tsv')
        for line in file_handle:
            gene_id, organism = line.rstrip().split("\t")
            false_gene_organism[gene_id] = organism

        false_gff = gff_file.HandleGFF(false_gff_file, false_gene_organism)
        false_gff.read_gff_file()

        false_gff.gene_organism = false_gene_organism

        false_gff.scan_gff_for_errors()
        false_gff.scan_mrna_sequence(fasta_file='./input_files/simple_false.fasta')

        sort_order_list = ('all', 'owner', 'organism_name', 'gene_id', 'mrna_id')

        error_object_list = list()
        error_lookup_table = dict()

        for key, value in false_gff.errors.items():
            error_object_list.append(value)

        error_lookup_table['all'] = error_object_list
        error_lookup_table['owner'] = list()
        error_lookup_table['organism_name'] = list()
        error_lookup_table['gene_id'] = list()
        error_lookup_table['mrna_id'] = list()

        annotation_quality_report.sort_and_write_errors(error_lookup_table, sort_order_list, 0, './')

        time_stamp = str(datetime.datetime.now().date())
        file_name = './' + 'simple@ebi.ac.uk' + '_' + time_stamp + '.error'
        self.assertEqual(True, filecmp.cmp('./output_files/simple_test_output.txt', file_name, shallow=False))

    def test_two_annotators(self):
        two_genes_false_gff = './input_files/two_annotator_false.gff'
        false_gene_organism = dict()
        file_handle = open('./input_files/two_gene_organism.tsv')
        for line in file_handle:
            gene_id, organism = line.rstrip().split("\t")
            false_gene_organism[gene_id] = organism

        false_gff = gff_file.HandleGFF(two_genes_false_gff, false_gene_organism)
        false_gff.read_gff_file()

        false_gff.gene_organism = false_gene_organism

        false_gff.scan_gff_for_errors()
        #false_gff.scan_mrna_sequence(fasta_file='./simple_false.fasta')

        sort_order_list = ('all', 'owner', 'organism_name', 'gene_id', 'mrna_id')

        error_object_list = list()
        error_lookup_table = dict()

        for key, value in false_gff.errors.items():
            error_object_list.append(value)
        error_lookup_table['all'] = error_object_list
        error_lookup_table['owner'] = list()
        error_lookup_table['organism_name'] = list()
        error_lookup_table['gene_id'] = list()
        error_lookup_table['mrna_id'] = list()

        annotation_quality_report.sort_and_write_errors(error_lookup_table, sort_order_list, 0, './')

        time_stamp = str(datetime.datetime.now().date())
        file_name1 = './' + 'annotator1@ebi.ac.uk' + '_' + time_stamp + '.error'
        self.assertEqual(True, filecmp.cmp('./output_files/annotator1_test_output.txt', file_name1, shallow=False))

        file_name2 = './' + 'annotator2@ebi.ac.uk' + '_' + time_stamp + '.error'
        self.assertEqual(True, filecmp.cmp('./output_files/annotator2_test_output.txt', file_name2, shallow=False))

    def test_two_species(self):

        two_genes_false_gff = './input_files/two_species_false.gff'
        false_gene_organism = dict()
        file_handle = open('./input_files/two_species_organism.tsv')
        for line in file_handle:
            gene_id, organism = line.rstrip().split("\t")
            false_gene_organism[gene_id] = organism

        false_gff = gff_file.HandleGFF(two_genes_false_gff, false_gene_organism)
        false_gff.read_gff_file()

        false_gff.gene_organism = false_gene_organism

        false_gff.scan_gff_for_errors()
        #false_gff.scan_mrna_sequence(fasta_file='./simple_false.fasta')

        sort_order_list = ('all', 'owner', 'organism_name', 'gene_id', 'mrna_id')

        error_object_list = list()
        error_lookup_table = dict()

        for key, value in false_gff.errors.items():
            error_object_list.append(value)
        error_lookup_table['all'] = error_object_list
        error_lookup_table['owner'] = list()
        error_lookup_table['organism_name'] = list()
        error_lookup_table['gene_id'] = list()
        error_lookup_table['mrna_id'] = list()

        annotation_quality_report.sort_and_write_errors(error_lookup_table, sort_order_list, 0, './')

        time_stamp = str(datetime.datetime.now().date())
        file_name = './' + 'species@ebi.ac.uk' + '_' + time_stamp + '.error'
        self.assertEqual(True, filecmp.cmp('./output_files/species_test_output.txt', file_name, shallow=False))

    def test_two_genes(self):

        two_genes_false_gff = './input_files/two_genes_false.gff'
        false_gene_organism = dict()
        file_handle = open('./input_files/two_gene_organism.tsv')
        for line in file_handle:
            gene_id, organism = line.rstrip().split("\t")
            false_gene_organism[gene_id] = organism

        false_gff = gff_file.HandleGFF(two_genes_false_gff, false_gene_organism)
        false_gff.read_gff_file()

        false_gff.gene_organism = false_gene_organism

        false_gff.scan_gff_for_errors()
        #false_gff.scan_mrna_sequence(fasta_file='./simple_false.fasta')

        sort_order_list = ('all', 'owner', 'organism_name', 'gene_id', 'mrna_id')

        error_object_list = list()
        error_lookup_table = dict()

        for key, value in false_gff.errors.items():
            error_object_list.append(value)
        error_lookup_table['all'] = error_object_list
        error_lookup_table['owner'] = list()
        error_lookup_table['organism_name'] = list()
        error_lookup_table['gene_id'] = list()
        error_lookup_table['mrna_id'] = list()

        annotation_quality_report.sort_and_write_errors(error_lookup_table, sort_order_list, 0, './')

        time_stamp = str(datetime.datetime.now().date())
        file_name = './' + 'twoGenes@ebi.ac.uk' + '_' + time_stamp + '.error'
        self.assertEqual(True, filecmp.cmp('./output_files/twoGenes_test_output.txt', file_name, shallow=False))

    def test_two_mrna(self):

        two_genes_false_gff = './input_files/two_mrna_false.gff'
        false_gene_organism = dict()
        file_handle = open('./input_files/two_mrna_organism.tsv')
        for line in file_handle:
            gene_id, organism = line.rstrip().split("\t")
            false_gene_organism[gene_id] = organism

        false_gff = gff_file.HandleGFF(two_genes_false_gff, false_gene_organism)
        false_gff.read_gff_file()

        false_gff.gene_organism = false_gene_organism

        false_gff.scan_gff_for_errors()
        #false_gff.scan_mrna_sequence(fasta_file='./simple_false.fasta')

        sort_order_list = ('all', 'owner', 'organism_name', 'gene_id', 'mrna_id')

        error_object_list = list()
        error_lookup_table = dict()

        for key, value in false_gff.errors.items():
            error_object_list.append(value)
        error_lookup_table['all'] = error_object_list
        error_lookup_table['owner'] = list()
        error_lookup_table['organism_name'] = list()
        error_lookup_table['gene_id'] = list()
        error_lookup_table['mrna_id'] = list()

        annotation_quality_report.sort_and_write_errors(error_lookup_table, sort_order_list, 0, './')

        time_stamp = str(datetime.datetime.now().date())
        file_name = './' + 'twomrna@ebi.ac.uk' + '_' + time_stamp + '.error'
        self.assertEqual(True, filecmp.cmp('./output_files/twomrna_test_output.txt', file_name, shallow=False))


if __name__ == '__main__':
    unittest.main()
