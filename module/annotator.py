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


class AnnotatorSummary:
    def __init__(self, email):
        self.email = email

        self.total_mrna_count = int()
        self.finished_mrna_count = int()

        self.total_gene_count = int()
        self.finished_gene_count = int()
        self.gene_list = list()
        self.unfinished_gene_list = list()

        self.total_pseudogene_count = int()
        self.finished_pseudogene_count = int()
        self.pseudogene_list = list()
        self.unfinished_pseudogene_list = list()

        self.non_canonical_count = int()

    def add_gene(self, name, finished=False):

        self.total_gene_count += 1
        if finished:
            self.finished_gene_count += 1
        else:
            self.unfinished_gene_list.append(name)
        self.gene_list.append(name)

    def add_pseudogene(self, name, finished=False):

        self.total_pseudogene_count += 1
        if finished:
            self.finished_pseudogene_count += 1
        else:
            self.unfinished_pseudogene_list.append(name)
        self.pseudogene_list.append(name)

    def add_mrna(self, finished=True):

        self.total_mrna_count += 1
        if finished:
            self.finished_mrna_count += 1

    def add_non_canonical(self):
        self.non_canonical_count += 1
