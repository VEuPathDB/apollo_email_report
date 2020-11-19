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

import gzip
import requests

from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional

from Bio import SeqIO # type: ignore
# NB: the above mentioned ^^^^^^^^ "type ignore" means "don't check this dep/package"

class GeneticCode:
    # this class can be moved somewhere else
    table: Dict[int, Dict[str, List[str]]] = {
            1 : {
                "START" : ["ATG"],
                "STOP" : ["TAA", "TAG", "TGA"]
            },
            # ...
        }

    @classmethod
    def starts(cls, code):
        if not self.table or code not in self.table:
            return []
        return self.table[code].get("START", [])

    @classmethod
    def stops(cls, code):
        if not self.table or code not in self.table:
            return []
        return self.table[code].get("STOP", [])


@dataclass
class CodingSequence:
    # may be usinf `dataclass` is an overkill in this case
    # internal, not added as init parameters
    sequence_type: str = field(init = False, default = 'cds')
    errors: Dict[str, str] = field(init = False, repr = False, default_factory = dict)
    valid_starts: FrozenSet[str] = field(init = False, default_factory = frozenset)
    valid_stops: FrozenSet[str] = field(init = False, default_factory = frozenset)
    # fields 
    feature_name: str
    organism_name: str
    sequence_name: str
    sequence: Optional[str] = field(repr = False, default = None)
    genetic_code: Optional[int] = field(default = 1)

    def __post_init__(self):
        self.valid_starts = fozenset(GeneticCode.starts(self.genetic_code))
        self.valid_stops = fozenset(GeneticCode.stops(self.genetic_code))

    def get_sequence_from_file(self, fasta_file: str) -> bool:
        try:
          _open = fasta_file.endswith(".gz") and gzip.open or open
          with _open(fasta_file, 'rt') as fasta:
            fasta_parser = SeqIO.parse(fasta, "fasta")
            for rec in fasta_parser:
              # m.b check if rec.name == self.feature_name
              self.sequence = rec.seq
              break
        except:
            return False

        if not self.sequence: return False
        return True

    def get_sequence_from_url(self, base_url: str) -> bool:
        try:
            url = base_url + "sequence/{}/{}/{}.{}?ignoreCache=true".format(
                self.organism_name,
                self.sequence_name,
                self.feature_name,
                CodingSequence.sequence_type
            )
            response = requests.get(url)
            if response and response.status_code == requests.codes.ok:
                self.sequence = response.text
        except:
            return False

        if not self.sequence: return False
        return True

    def get_sequence(self, uri:str) -> bool:
        return self.get_sequence_from_file(uri) or self.get_sequence_from_url(uri)

    def error(self, tag: str, message: str) -> bool:
        self.errors[tag] = message
        return False

    def has_start_codon(self):
        seq = self.sequence
        if not seq or len(seq) < 3: return False
        first_codon = seq[0:3]
        if first_codon in self.valid_starts:
            return True
        return self.error('start_codon', 'no start codon')

    def has_stop_codon(self):
        seq = self.sequence
        if not seq or len(seq) < 3: return False
        last_codon = self.sequence[-3:]
        if last_codon in self.valid_stops:
            return True
        return self.error('stop_codon', 'no stop codon')

    def no_internal_stop_codon(self):
        seq = self.sequence
        if not seq or len(seq) < 3: return True
        internal_stops = 0
        for i in range(0, len(seq) - 3, 3):
            codon = seq[i:i + 3]
            if codon in self.valid_stops:
                internal_stops += 1
        if internal_stops:
            return self.error('no_internal_stop_codon', '%d internal stop codon' % (internal_stops))
        return True

