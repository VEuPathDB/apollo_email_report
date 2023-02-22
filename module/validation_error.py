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


class ValidationError:
    def __init__(
        self,
        owner: str,
        organism_name: str,
        gene_id: str,
        mrna_id: str = None,
        gene_name: str = None,
        locus: str = None,
    ) -> None:
        self.owner = owner
        self.organism_name = organism_name
        self.gene_id = gene_id
        self.mrna_id = mrna_id
        self.gene_name = gene_name
        self.locus = locus
        self.gff_format_error = list()
        self.sequence_error = list()

    def add_gff_format_error(
        self,
        field_type: str,
        feature_type: str,
        feature_id: str,
        feature_value: str,
        parent_id: str,
        parent_value: str,
    ) -> None:
        error_values = {
            "field_type": field_type,
            "feature_type": feature_type,
            "feature_id": feature_id,
            "feature_value": feature_value,
            "parent_id": parent_id,
            "parent_value": parent_value,
        }
        self.gff_format_error.append(error_values)

    def add_sequence_error(
        self, mrna_id: str, error_name: str, error_text: str
    ) -> None:
        error_values = {
            "mrna_id": mrna_id,
            "error_name": error_name,
            "error_text": error_text,
        }
        self.sequence_error.append(error_values)

    def gff_error_text(self) -> List:
        output_list = list()
        for error_values in self.gff_format_error:
            output_text = (
                "The {field_type} of this feature {feature_type} : {feature_id} with value"
                " {feature_value} is not in accordance with"
                " its parent {parent_id} with value {parent_value}".format(
                    **error_values
                )
            )
            output_list.append(output_text)
        return output_list

    def sequence_error_text(self) -> List:
        output_list = list()
        for error_values in self.sequence_error:
            output_text = (
                "The coding sequence for mRNA: {mrna_id} has {error_text}".format(
                    **error_values
                )
            )
            output_list.append(output_text)
        return output_list
