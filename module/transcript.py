"""
Copyright [2017-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import requests
import re
import urllib.parse


class CodingSequence:
    sequence_type = 'cds'

    def __init__(self, feature_name, organism_name, sequence_name):
        self.feature_name = feature_name
        self.organism_name = organism_name
        self.sequence_name = sequence_name
        self.errors = dict()
        self.sequence = str()

    def get_sequence(self, base_url=None, username=None, password=None, fasta_file=None):

        if base_url:
            url = urllib.parse.urljoin(base_url, 'sequence/sequenceByName')
            body = {'username': username, 'password': password, 'organismString': self.organism_name,
                    'sequenceName': self.sequence_name, 'featureName': self.feature_name,
                    'type': CodingSequence.sequence_type, 'ignoreCache':'true'}

            response = requests.post(url, json=body)
            seq = response.text
            if response.status_code == requests.codes.ok:
                matches = re.match(r'^([CTGAN]+)<', seq)
                if matches:
                    print(f"Remove weird html part at the end of {self.feature_name}")
                    print(f"[{seq}]")
                    seq = matches[1]
                if not re.match(r'^[CGTAN]+$', seq):
                    print(f"Incorrect CDS sequence: [{seq}]")
                    return False
                self.sequence = seq
            else:
                print(f"No sequence retrieved for {self.feature_name}: code {response.status_code}")
                return False
        elif fasta_file:
            with open(fasta_file, 'r') as file_handle:
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
