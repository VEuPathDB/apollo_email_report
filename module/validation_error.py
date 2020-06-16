class ValidationError:
    def __init__(self, owner, organism_name, gene_id, mrna_id=None, gene_name=None, locus=None):
        self.owner = owner
        self.organism_name = organism_name
        self.gene_id = gene_id
        self.mrna_id = mrna_id
        self.gene_name = gene_name
        self.locus = locus
        self.gff_format_error = list()
        self.sequence_error = list()

    def add_gff_format_error(self, field_type, feature_type, feature_id, feature_value, parent_id, parent_value):
        error_values = {'field_type': field_type, 'feature_type': feature_type, 'feature_id': feature_id,
                        'feature_value': feature_value, 'parent_id': parent_id, 'parent_value': parent_value}
        self.gff_format_error.append(error_values)

    def add_sequence_error(self, mrna_id, error_name, error_text):
        error_values = {'mrna_id': mrna_id, 'error_name': error_name, 'error_text': error_text}
        self.sequence_error.append(error_values)

    def gff_error_text(self):
        output_list = list()
        for error_values in self.gff_format_error:

            output_text = "The {field_type} of this feature {feature_type} : {feature_id} with value" \
                          " {feature_value} is not in accordance with" \
                          " its parent {parent_id} with value {parent_value}.\n".format(**error_values)
            output_list.append(output_text)
        return output_list

    def sequence_error_text(self):
        output_list = list()
        for error_values in self.sequence_error:
            output_text = "The coding sequence for mRNA: {mrna_id} has {error_text}.\n".format(**error_values)
            output_list.append(output_text)
        return output_list
