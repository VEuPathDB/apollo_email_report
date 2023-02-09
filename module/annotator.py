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
        self.mrnas = list()
        self.unfinished_mrnas = list()

        self.total_gene_count = int()
        self.finished_gene_count = int()
        self.genes = list()
        self.unfinished_genes = list()

        self.total_pseudogene_count = int()
        self.finished_pseudogene_count = int()
        self.pseudogenes = list()
        self.unfinished_pseudogenes = list()

        self.total_ncrna_count = int()
        self.finished_ncrna_count = int()
        self.ncrnas = list()
        self.unfinished_ncrnas = list()

        self.non_canonical_count = int()

    def add_gene(self, name, finished=False):

        self.total_gene_count += 1
        if finished:
            self.finished_gene_count += 1
        else:
            self.unfinished_genes.append(name)
        self.genes.append(name)

    def add_pseudogene(self, name, finished=False):

        self.total_pseudogene_count += 1
        if finished:
            self.finished_pseudogene_count += 1
        else:
            self.unfinished_pseudogenes.append(name)
        self.pseudogenes.append(name)

    def add_ncrna(self, name, finished=False):

        self.total_ncrna_count += 1
        if finished:
            self.finished_ncrna_count += 1
        else:
            self.unfinished_ncrnas.append(name)
        self.ncrnas.append(name)

    def add_mrna(self, name, finished=True):

        self.total_mrna_count += 1
        if finished:
            self.finished_mrna_count += 1
        else:
            self.unfinished_mrnas.append(name)
        self.mrnas.append(name)

    def add_non_canonical(self):
        self.non_canonical_count += 1

    def has_changes(self) -> bool:
        return (self.total_gene_count + self.total_pseudogene_count > 0)

    def has_unfinished(self) -> bool:
        total = self.total_gene_count + self.total_pseudogene_count
        finished = self.finished_gene_count + self.finished_pseudogene_count
        return total - finished > 0

    def get_unfinished(self) -> List[str]:
        unfinished = []
        unfinished += [f"gene\t{name}" for name in self.unfinished_genes]
        unfinished += [f"mRNA\t{name}" for name in self.unfinished_mrnas]
        unfinished += [f"ncRNA\t{name}" for name in self.unfinished_ncrnas]
        unfinished += [f"pseudogene\t{name}" for name in self.unfinished_pseudogenes]
        return unfinished

    def get_all(self) -> List[str]:
        return self.genes