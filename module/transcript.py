import requests


class CodingSequence:
    sequence_type = 'cds'

    def __init__(self, feature_name, organism_name, sequence_name):
        self.feature_name = feature_name
        self.organism_name = organism_name
        self.sequence_name = sequence_name
        self.errors = dict()
        self.sequence = str()

    def get_sequence(self, base_url=None, fasta_file=None):

        if base_url:
            url = base_url + "sequence/{}/{}/{}.{}?ignoreCache=true".format(self.organism_name, self.sequence_name,
                                                                            self.feature_name,
                                                                            CodingSequence.sequence_type)
            response = requests.get(url)
            seq = response.text
            if response.status_code == requests.codes.ok:
                self.sequence = seq
            else:
                return False
        elif fasta_file:
            file_handle = open(fasta_file, 'r')
            for seq in file_handle:
                if seq[0] != '>':
                    self.sequence += seq.rstrip()
        else:
            return False

    def coding_sequence_has_start_codon(self):
        if self.sequence[0:3] == 'ATG':
            return True
        else:
            self.errors['start_codon'] = 'no start codon'
            return False

    def coding_sequence_has_stop_codon(self):
        last_codon = self.sequence[-3:]
        if last_codon == 'TAA' or last_codon == 'TAG' or last_codon == 'TGA':
            return True
        else:
            self.errors['stop_codon'] = 'no stop codon'
            return False

    def coding_sequence_no_internal_stop_codon(self):
        internal_stop_codon_count = 0
        for i in range(0, len(self.sequence) - 3, 3):
            codon = self.sequence[i:i + 3]
            if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                internal_stop_codon_count += 1
                self.errors['no_internal_stop_codon'] = str(internal_stop_codon_count) + ' internal stop codon'
        if internal_stop_codon_count:
            return False
        else:
            return True
