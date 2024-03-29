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


from typing import List


class AnnotatorSummary:
    def __init__(self, email):
        self.email = email

        self.total_mrna_count = int()
        self.finished_mrna_count = int()
        self.mrnas = set()
        self.unfinished_mrnas = set()

        self.total_gene_count = int()
        self.finished_gene_count = int()
        self.genes = set()
        self.unfinished_genes = set()

        self.total_pseudogene_count = int()
        self.finished_pseudogene_count = int()
        self.pseudogenes = set()
        self.unfinished_pseudogenes = set()

        self.total_ncrna_count = int()
        self.finished_ncrna_count = int()
        self.ncrnas = set()
        self.unfinished_ncrnas = set()

        self.non_canonical_count = int()

    def add_gene(self, name, finished=False):
        self.genes.add(name)
        if not finished:
            self.unfinished_genes.add(name)
        self.total_gene_count = len(self.genes)
        self.finished_gene_count = self.total_gene_count - len(self.unfinished_genes)

    def add_pseudogene(self, name, finished=False):
        self.pseudogenes.add(name)
        if not finished:
            self.unfinished_pseudogenes.add(name)
        self.total_pseudogene_count = len(self.pseudogenes)
        self.finished_pseudogene_count = self.total_pseudogene_count - len(
            self.unfinished_pseudogenes
        )

    def add_ncrna(self, name, finished=False):
        self.ncrnas.add(name)
        if not finished:
            self.unfinished_ncrnas.add(name)
        self.total_ncrna_count = len(self.ncrnas)
        self.finished_ncrna_count = self.total_ncrna_count - len(self.unfinished_ncrnas)

    def add_mrna(self, name, finished=True):
        self.mrnas.add(name)
        if not finished:
            self.unfinished_mrnas.add(name)
        self.total_mrna_count = len(self.mrnas)
        self.finished_mrna_count = self.total_mrna_count - len(self.unfinished_mrnas)

    def add_non_canonical(self):
        self.non_canonical_count += 1

    def has_changes(self) -> bool:
        return self.total_gene_count + self.total_pseudogene_count > 0

    def has_unfinished(self) -> bool:
        total = self.total_gene_count + self.total_pseudogene_count
        finished = self.finished_gene_count + self.finished_pseudogene_count
        return total - finished > 0

    def get_unfinished(self) -> List[str]:
        unfinished = []
        unfinished += [
            f"protein_coding\t{name}" for name in sorted(self.unfinished_mrnas)
        ]
        unfinished += [f"ncRNA\t{name}" for name in sorted(self.unfinished_ncrnas)]
        unfinished += [
            f"pseudogene\t{name}" for name in sorted(self.unfinished_pseudogenes)
        ]
        return unfinished
