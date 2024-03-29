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
from typing import Dict, List
from module import annotator, transcript, validation_error

any_organism = "ANY_ORGANISM"
MAX_ITERATION = 6

top_level_feat = set(("gene", "pseudogene"))
sub_prot_feat = set(
    (
        "mRNA",
        "exon",
        "CDS",
        "non_canonical_five_prime_splice_site",
        "non_canonical_three_prime_splice_site",
    )
)
non_coding_feat = set(
    (
        "rRNA",
        "tRNA",
        "sRNA",
        "snoRNA",
        "ncRNA",
        "miRNA",
        "guide_RNA",
        "RNAse_MRP_RNA",
        "telomerase_RNA",
        "SRP_RNA",
        "lnc_RNA",
        "RNAse_P_RNA",
        "scRNA",
        "piRNA",
        "tmRNA",
        "enzymatic_RNA",
    )
)
allowed_feature = top_level_feat.union(sub_prot_feat).union(non_coding_feat)


class HandleGFF:
    def __init__(self, file_path, gene_organism, moderator):
        self.file_path = file_path
        self.gene_organism = gene_organism
        self.moderator = moderator
        self.feature_type = dict()
        self.gene_meta_info = dict()
        self.feature_owner = dict()
        self.child_parent_relationship = dict()
        self.transcripts = list()
        self.annotators = dict()
        self.parent_name = dict()

        self.fields = dict()
        self.errors = dict()

    def read_gff_file(self):
        with open(self.file_path, "r") as file_handle:
            line_number = 0
            disqualified_features = list()
            finished_gene_status = False
            for line in file_handle:
                fields = line.rstrip().split("\t")
                line_number += 1
                if len(fields) != 9:
                    finished_gene_status = False
                    continue  # skip line as not GFF
                (
                    feature_type,
                    owner,
                    scaffold,
                    strand,
                    feature_id,
                    parent_id,
                    name,
                    locus,
                    status,
                    partial,
                ) = extract_fields_from_gff(fields)
                if owner is not None:
                    owner = self.get_first_owner(owner)
                if feature_type not in allowed_feature:
                    disqualified_features.append(feature_id)

                if feature_id in disqualified_features:
                    continue
                if parent_id in disqualified_features:
                    disqualified_features.append(feature_id)
                    continue

                if feature_type == "gene":
                    finished_gene_status = False
                    self.gene_meta_info[feature_id] = (name, locus)
                    if owner is not None:
                        if owner not in self.annotators:
                            self.annotators[owner] = annotator.AnnotatorSummary(owner)

                        if status == "Finished" or status == "Finished annotating":
                            finished_gene_status = True
                        self.annotators[owner].add_gene(name, finished_gene_status)
                    else:
                        owner = self.moderator
                        print("No owner for Gene: " + feature_id)

                elif feature_type == "pseudogene":
                    finished_gene_status = False
                    self.gene_meta_info[feature_id] = (name, locus)
                    if owner is not None:
                        if owner not in self.annotators:
                            self.annotators[owner] = annotator.AnnotatorSummary(owner)
                        if status == "Finished" or status == "Finished annotating":
                            finished_gene_status = True
                        self.annotators[owner].add_pseudogene(
                            name, finished_gene_status
                        )
                    else:
                        owner = self.moderator
                        print("No owner for Pseudogene: " + feature_id)

                elif feature_type == "mRNA":
                    organism = self.gene_organism.get(parent_id)
                    if not organism:
                        organism = any_organism

                    if owner is not None:
                        if owner not in self.annotators:
                            self.annotators[owner] = annotator.AnnotatorSummary(owner)

                        parent_name, _ = self.gene_meta_info[parent_id]
                        if finished_gene_status:
                            self.annotators[owner].add_mrna(parent_name, True)

                            # Save the transcript for checking
                            if not partial:
                                self.transcripts.append(
                                    (feature_id, organism, scaffold)
                                )
                        else:
                            self.annotators[owner].add_mrna(parent_name, False)
                    else:
                        owner = self.moderator
                        print("No owner for mRNA: " + feature_id)

                elif feature_type in non_coding_feat:
                    organism = self.gene_organism.get(parent_id)
                    if not organism:
                        organism = any_organism

                    if owner is not None:
                        if owner not in self.annotators:
                            self.annotators[owner] = annotator.AnnotatorSummary(owner)

                        parent_name, _ = self.gene_meta_info[parent_id]
                        self.annotators[owner].add_ncrna(
                            parent_name, finished_gene_status
                        )
                    else:
                        owner = self.moderator
                        print("No owner for ncRNA: " + feature_id)

                if feature_id in self.child_parent_relationship:
                    feature_id = feature_id + "_" + str(line_number)

                if feature_type in [
                    "non_canonical_five_prime_splice_site",
                    "non_canonical_three_prime_splice_site",
                ]:
                    owner = self.get_parent_owner(parent_id)
                    self.annotators[owner].add_non_canonical()

                self.feature_owner[feature_id] = owner
                self.feature_type[feature_id] = feature_type

                if finished_gene_status:
                    self.fields[("scaffold", feature_id)] = scaffold
                    self.fields[("strand", feature_id)] = strand
                    self.fields[("position", feature_id)] = (
                        int(fields[3]),
                        int(fields[4]),
                    )

                if parent_id:
                    self.child_parent_relationship[feature_id] = parent_id
                elif feature_type in top_level_feat:
                    self.child_parent_relationship[
                        feature_id
                    ] = feature_id  # to avoid checking if ID exists
                else:
                    print("Feature not recognized", feature_type, feature_id)

    @staticmethod
    def get_first_owner(owner):
        owners = owner.rstrip().split(",")
        return owners[0]

    def get_gene_id(self, feature_id, iteration=0):
        if iteration > MAX_ITERATION:
            raise RecursionError(f"Too many iteration for id {feature_id}")
        if feature_id == self.child_parent_relationship[feature_id]:
            return feature_id
        else:
            parent_id = self.child_parent_relationship[feature_id]
            iteration += 1
            return self.get_gene_id(parent_id, iteration)

    def get_parent_owner(self, feature_id, iteration=0):
        if iteration > MAX_ITERATION:
            raise RecursionError(f"Too many iteration for id {feature_id}")
        if self.feature_owner[feature_id]:
            return self.feature_owner[feature_id]
        else:
            parent_id = self.child_parent_relationship[feature_id]
            iteration += 1
            return self.get_parent_owner(parent_id)

    def get_mrna_id(self, feature_id, iteration=0):
        if iteration > MAX_ITERATION:
            raise RecursionError(f"Too many iteration for id {feature_id}")
        if self.feature_type[feature_id] == "mRNA":
            return feature_id
        elif self.feature_type[feature_id] in top_level_feat:
            return None
        else:
            parent_id = self.child_parent_relationship[feature_id]
            iteration += 1
            return self.get_mrna_id(parent_id, iteration)

    def scan_gff_for_errors(self):
        for key, value in self.fields.items():
            field, feature_id = key

            if feature_id not in self.child_parent_relationship:
                continue

            parent_id = self.child_parent_relationship[feature_id]
            parent_value = self.fields[(field, parent_id)]

            gene_id = self.get_gene_id(feature_id)
            owner = self.get_parent_owner(feature_id)
            mrna_id = self.get_mrna_id(feature_id)

            gene_name, locus = self.gene_meta_info[gene_id]
            arguments = dict()
            arguments["owner"] = owner
            arguments["organism_name"] = self.gene_organism.get(gene_id)
            if not arguments["organism_name"]:
                arguments["organism_name"] = any_organism
            arguments["gene_id"] = gene_id
            arguments["mrna_id"] = mrna_id
            arguments["gene_name"] = gene_name
            arguments["locus"] = locus

            if field == "position":
                begin, end = value
                parent_begin, parent_end = parent_value
                if not (
                    parent_begin <= begin <= parent_end
                    and parent_begin <= end <= parent_end
                ):
                    arguments["field_type"] = field
                    arguments["feature_type"] = self.feature_type[feature_id]
                    arguments["feature_value"] = "Begin:{}..End:{}".format(begin, end)
                    arguments["parent_id"] = parent_id
                    arguments["parent_value"] = "Begin:{}..End:{}".format(
                        parent_begin, parent_end
                    )

                    self.add_validation_error(feature_id, "gff_error", **arguments)

            else:
                if value != parent_value:
                    arguments["field_type"] = field
                    arguments["feature_type"] = self.feature_type[feature_id]
                    arguments["feature_value"] = value
                    arguments["parent_id"] = parent_id
                    arguments["parent_value"] = parent_value
                    self.add_validation_error(feature_id, "gff_error", **arguments)

    def scan_mrna_sequence(
        self, base_url=None, username=None, password=None, fasta_file=None
    ):
        for mrna_id, organism, scaffold in self.transcripts:
            mrna = transcript.CodingSequence(mrna_id, organism, scaffold)

            if base_url:
                mrna.get_sequence(
                    base_url=base_url, username=username, password=password
                )
            elif fasta_file:
                mrna.get_sequence(fasta_file=fasta_file)
            else:
                return False

            if not mrna.sequence:
                print(f"No mRNA sequence could be used for {mrna_id}")
                continue

            mrna.coding_sequence_has_start_codon()
            mrna.coding_sequence_has_stop_codon()
            mrna.coding_sequence_no_internal_stop_codon()

            if mrna.errors != {}:
                gene_id = self.get_gene_id(mrna_id)
                owner = self.get_parent_owner(mrna_id)

                gene_name, locus = self.gene_meta_info[gene_id]
                arguments = dict()
                arguments["owner"] = owner
                arguments["organism_name"] = self.gene_organism.get(gene_id)
                if not arguments["organism_name"]:
                    arguments["organism_name"] = any_organism
                arguments["gene_id"] = gene_id
                arguments["mrna_id"] = mrna_id
                arguments["gene_name"] = gene_name
                arguments["locus"] = locus

                for key, value in mrna.errors.items():
                    arguments["error_name"] = key
                    arguments["error_text"] = value
                    self.add_validation_error(mrna_id, "sequence_error", **arguments)

    def add_validation_error(self, feature_id, error_type, **kwargs):
        error_id = feature_id
        if error_id not in self.errors:
            self.errors[error_id] = validation_error.ValidationError(
                kwargs["owner"],
                kwargs["organism_name"],
                kwargs["gene_id"],
                kwargs["mrna_id"],
                kwargs["gene_name"],
                kwargs["locus"],
            )

        if error_type == "sequence_error":
            self.errors[error_id].add_sequence_error(
                kwargs["mrna_id"], kwargs["error_name"], kwargs["error_text"]
            )
        elif error_type == "gff_error":
            self.errors[error_id].add_gff_format_error(
                kwargs["field_type"],
                kwargs["feature_type"],
                feature_id,
                kwargs["feature_value"],
                kwargs["parent_id"],
                kwargs["parent_value"],
            )


def parse_attribs(attribs_field: str) -> Dict:
    """Returns the attribs as a Dict."""
    attribs_str = attribs_field.split(";")

    attribs = {}
    for pair in attribs_str:
        (key, value) = pair.split("=")
        attribs[key] = value

    return attribs


def extract_fields_from_gff(fields: List) -> Dict:
    """Return a list of fields from a gff feature line"""
    scaffold, _, feature_type, begin, end, _, strand, _, attribs_field = fields
    attribs: dict = parse_attribs(attribs_field)

    feature_id = attribs.get("ID")
    parent_id = attribs.get("Parent")
    owner = attribs.get("owner")
    name = attribs.get("Name")
    status = attribs.get("status")
    partial = False
    if (
        attribs.get("is_fmax_partial", "") == "true"
        or attribs.get("is_fmin_partial", "") == "true"
    ):
        partial = True

    locus = None
    if feature_type == "gene":
        locus = "{}:{}..{}".format(scaffold, begin, end)

    if feature_id:
        return (
            feature_type,
            owner,
            scaffold,
            strand,
            feature_id,
            parent_id,
            name,
            locus,
            status,
            partial,
        )
    else:
        return False
