class AnnotatorSummary:
    def __init__(self, email):
        self.email = email
        self.total_gene_count = int()
        self.total_mrna_count = int()
        self.finished_gene_count = int()
        self.finished_mrna_count = int()
        self.gene_list = list()
        self.unfinished_gene_list = list()

    def add_gene(self, name, finished=False):

        self.total_gene_count += 1
        if finished:
            self.finished_gene_count += 1
        else:
            self.unfinished_gene_list.append(name)
        self.gene_list.append(name)

    def add_mrna(self, finished=False):

        self.total_mrna_count += 1
        if finished:
            self.finished_mrna_count += 1
