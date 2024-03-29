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
import inspect
import json
import os
import sys
from pathlib import Path
import shutil

from module import annotation_quality_report as report, gff_file
from module.gff_file import HandleGFF

root_path = Path(inspect.getfile(sys.modules[__name__])).parent / ".."
summary_footer_path = root_path / "summary_static_footer.txt"
error_footer_path = root_path / "error_static_footer.txt"


class ApolloReporter:
    def __init__(self, config) -> None:
        self.config = config

    def prepare_summary_gff(self):
        config = self.config
        gff_dir = config["SETUP"]["dir"]
        organism_file = config["SETUP"]["organism_file"]
        self._make_new_directory(gff_dir)
        organism_file_handle = open(organism_file)
        master_gff_file_name = gff_dir + "/master.gff"
        master_gff_file_handle = open(master_gff_file_name, "w")

        gene_to_organism_lookup = dict()
        for organism_line in organism_file_handle:
            organism = organism_line.rstrip()
            gff_file_path = self.download_organism_gff(organism, gff_dir)

            if not gff_file_path:
                continue

            self._join_gff_to_master(
                organism, gff_file_path, master_gff_file_handle, gene_to_organism_lookup
            )

        master_gff_file_handle.close()
        organism_lookup_file_name = gff_dir + "/organism_lookup.json"
        with open(organism_lookup_file_name, "w") as outfile:
            json.dump(gene_to_organism_lookup, outfile)

        return gene_to_organism_lookup, master_gff_file_name

    def download_organism_gff(self, organism_name, out_dir):
        config = self.config
        apollo_url = config["APOLLO"]["base_url"]
        apollo_user_name = config["APOLLO"]["username"]
        apollo_password = config["APOLLO"]["password"]
        clean_organism_name = organism_name.replace("/", "")
        organism_new_dir = out_dir + "/" + clean_organism_name + "/"
        self._make_new_directory(organism_new_dir)
        gff_file_name = report.download_gff(
            apollo_url,
            apollo_user_name,
            apollo_password,
            organism_name,
            organism_new_dir,
        )
        return gff_file_name

    def _join_gff_to_master(
        self, organism, gff_file_name, master_gff_file_handle, gene_organism
    ):
        gff_file_handle = open(gff_file_name, "r")
        for line in gff_file_handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
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
            ) = gff_file.extract_fields_from_gff(fields)

            if feature_type == "gene":
                gene_organism[feature_id] = organism

            master_gff_file_handle.write(line)

    @staticmethod
    def _make_new_directory(new_dir):
        if os.path.exists(new_dir):
            shutil.rmtree(new_dir)
            os.mkdir(new_dir)
        else:
            os.mkdir(new_dir)

    def prepare_recent_gff(self):
        config = self.config
        base_url = config["APOLLO"]["base_url"]
        username = config["APOLLO"]["username"]
        password = config["APOLLO"]["password"]
        days = int(config["SETUP"]["days"])
        out_dir = config["SETUP"]["dir"]

        self._make_new_directory(out_dir)

        recent_apollo_genes = report.get_recent_genes_from_apollo(
            base_url, username, password, days
        )
        gff_file_path = False
        if recent_apollo_genes:
            gff_file_path = report.get_gff(
                base_url, username, password, recent_apollo_genes, out_dir
            )
        else:
            print("No genes have been changed")
            exit()
        return recent_apollo_genes, gff_file_path

    def prepare_summary_emails(
        self, master_gff_file_name, gene_organism, file_extension
    ):
        config = self.config
        email_dir = Path(config["SETUP"]["dir"])
        footer_text = self.load_summary_footer()

        gff_file_object = gff_file.HandleGFF(
            master_gff_file_name, gene_organism, config["EMAIL"]["moderator"]
        )
        gff_file_object.read_gff_file()
        self._write_email_body(gff_file_object, email_dir)

        messages = list()
        annotation_summary_emails = self._collect_files(email_dir, file_extension)
        for email in annotation_summary_emails:
            email_address, email_message = self._compose_message(
                email_dir, email, footer_text
            )
            messages.append((email_address, email_message))
        return messages

    @staticmethod
    def _write_email_body(file_object, out_dir: Path):
        for owner, annotator_object in file_object.annotators.items():
            report.write_summary_text(annotator_object, out_dir)

    @staticmethod
    def _collect_files(out_dir: Path, filter_term):
        file_filtered_list = list()
        _, _, file_list = next(os.walk(str(out_dir)), (None, None, []))
        for file_name in file_list:
            _, file_extension = os.path.splitext(file_name)
            if "." + filter_term == file_extension:
                file_filtered_list.append(file_name)
        return file_filtered_list

    @staticmethod
    def _compose_message(out_dir: Path, email_body, footer):
        email_path = out_dir / email_body
        with email_path.open("r") as email_fh:
            address = email_fh.readline().rstrip()
            message = email_fh.readlines()
            message += footer

        return address, message

    def load_error_footer(self):
        config = self.config
        static_footer = config["EMAIL"].get("error_static_footer")
        if not static_footer:
            static_footer = error_footer_path
        else:
            static_footer = Path(static_footer)
        with static_footer.open("r") as footer_fh:
            return footer_fh.readlines()

    def load_summary_footer(self):
        config = self.config
        static_footer = config["EMAIL"].get("summary_static_footer")
        if not static_footer:
            static_footer = summary_footer_path
        else:
            static_footer = Path(static_footer)
        with static_footer.open("r") as footer_fh:
            return footer_fh.readlines()

    def prepare_error_emails(
        self, master_gff_file_name, gene_organism_lookup, file_extension
    ):
        config = self.config
        email_dir = Path(config["SETUP"]["dir"])
        apollo_url = config["APOLLO"]["base_url"]
        footer_text = self.load_error_footer()
        messages = list()

        gff_file_object = report.validate_gff(
            apollo_url,
            config["APOLLO"]["username"],
            config["APOLLO"]["password"],
            master_gff_file_name,
            gene_organism_lookup,
            config["EMAIL"]["moderator"],
        )

        if not gff_file_object:
            print("No error was found")
            return messages

        report.write_email_texts(gff_file_object.errors, email_dir)

        annotation_error_emails = self._collect_files(email_dir, file_extension)
        for email in annotation_error_emails:
            email_address, email_message = self._compose_message(
                email_dir, email, footer_text
            )
            messages.append((email_address, email_message))
        return messages

    @staticmethod
    def _sort_error_and_write_email_body(gff_file_object: HandleGFF, email_dir: Path) -> None:

        sort_order_list = ("all", "owner", "organism_name", "gene_id", "mrna_id")
        error_object_list = list()
        error_lookup_table = dict()

        for key, value in gff_file_object.errors.items():
            error_object_list.append(value)

        error_lookup_table["all"] = error_object_list
        error_lookup_table["owner"] = list()
        error_lookup_table["organism_name"] = list()
        error_lookup_table["gene_id"] = list()
        error_lookup_table["mrna_id"] = list()

        report.sort_and_write_errors_old(error_lookup_table, sort_order_list, 0, email_dir)

    def send_emails(self, email_type, list_of_emails):
        config = self.config
        email_dir = config["SETUP"]["dir"]
        mailgun_url = config["MAILGUN"]["url"]
        mailgun_key = config["MAILGUN"]["api_key"]
        from_address = config["EMAIL"]["from_address"]
        email_url = config["EMAIL"]["base_url"]
        client_id = config["EMAIL"]["client_id"]
        client_secret = config["EMAIL"]["client_secret"]

        mode = config["SETUP"]["mode"]

        if email_type == "summary":
            subject = config["EMAIL"]["summary_subject"]
        else:
            subject = config["EMAIL"]["error_subject"]

        for email in list_of_emails:
            user_id, email_message = email
            email_address = report.get_email(
                email_url, client_id, client_secret, user_id
            )
            moderators = config["EMAIL"]["moderator"]
            moderator_email_address = moderators
            if not email_address:
                email_address = config["EMAIL"]["moderator"]

            if email_type == "summary":
                file_attached = email_dir + "/" + user_id + ".gene_list"
            else:
                file_attached = None

            if mode != "live":
                email_address = config["EMAIL"]["moderator"]

            report.send_email_mailgun(
                mailgun_url,
                mailgun_key,
                from_address,
                email_address,
                moderator_email_address,
                subject,
                email_message,
                file_attached,
            )
