import unittest
from module import validation_error


class MyTestCase(unittest.TestCase):
    def test_GFF_format(self):
        expected_error_txt = (
            "The scaffold of this feature mrna : 222-xxx-111 with value"
            " chr1 is not in accordance with"
            " its parent 111-xxx-111 with value chr2.\n"
        )
        error = validation_error.ValidationError(
            "mikkel", "aedes", "111-xxx-111", "222-xxx-111"
        )
        error.add_gff_format_error(
            "scaffold", "mrna", "222-xxx-111", "chr1", "111-xxx-111", "chr2"
        )
        error_list = error.gff_error_text()

        for error_text in error_list:
            self.assertEqual(error_text, expected_error_txt)


if __name__ == "__main__":
    unittest.main()
