import unittest
from module import gff_file


class MyTestCase(unittest.TestCase):
    def test_summary(self):
        mini_gff_file = "./input_files/mini.gff"
        gene_organism = {
            "9519cfca-0c42-44d4-ab09-d37d33245d07": "sand_box",
            "78fbaf52-0a8b-4acb-a93c-5d50238e47bf": "sand_box",
        }

        mini_gff = gff_file.HandleGFF(mini_gff_file, gene_organism, "")
        mini_gff.read_gff_file()
        print(mini_gff.annotators)

        self.assertEqual(2, len(mini_gff.annotators))
        self.assertEqual(True, "annotator1@ebi.ac.uk" in mini_gff.annotators)
        self.assertEqual(True, "annotator2@ebi.ac.uk" in mini_gff.annotators)
        self.assertEqual(
            1, mini_gff.annotators["annotator1@ebi.ac.uk"].total_gene_count
        )
        self.assertEqual(
            1, mini_gff.annotators["annotator1@ebi.ac.uk"].finished_gene_count
        )
        self.assertEqual(
            2, mini_gff.annotators["annotator2@ebi.ac.uk"].finished_mrna_count
        )
        self.assertEqual(1, len(mini_gff.annotators["annotator2@ebi.ac.uk"].gene_list))


if __name__ == "__main__":
    unittest.main()
