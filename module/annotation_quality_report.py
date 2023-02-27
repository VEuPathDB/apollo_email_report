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
from pathlib import Path
from typing import Dict, List
import requests
from requests.exceptions import ConnectionError
import datetime
import urllib.parse
import re
from module import gff_file
from module.annotator import AnnotatorSummary
from module.validation_error import ValidationError


def get_recent_genes_from_apollo(base_url, username, password, days=1):
    webservice_data = {"username": username, "password": password, "days": days}
    url = urllib.parse.urljoin(base_url, "annotationEditor/getRecentAnnotations")

    response = requests.post(url, json=webservice_data)
    if response.status_code == requests.codes.ok:
        try:
            return response.json()
        except Exception:
            print(response.text)
            raise Exception()
    else:
        if response.status_code == requests.codes.unauthorized:
            raise Exception("Unauthorized")
        else:
            raise Exception(f"Response error: {response.status_code}")


def get_organisms(base_url, username, password):
    webservice_data = {"username": username, "password": password}
    url = urllib.parse.urljoin(base_url, "organism/findAllOrganisms")

    response = requests.post(url, json=webservice_data)
    if response.status_code == requests.codes.ok:
        try:
            return response.json()
        except Exception:
            print(response.text)
            raise Exception()
    else:
        if response.status_code == requests.codes.unauthorized:
            raise Exception("Unauthorized")
        else:
            raise Exception(f"Response error: {response.status_code}")


def get_email(base_url, client_id, client_secret, user_id):
    user_obj = re.match(r".*\.(\d+?)$", user_id, flags=0)

    if user_obj:
        upenn_id = int(user_obj.group(1))
    else:
        return False

    body = {
        "client_id": client_id,
        "client_secret": client_secret,
        "query": {"userId": upenn_id},
    }
    response = requests.post(base_url, json=body)
    if response.status_code == requests.codes.ok:
        return response.json()["email"]
    else:
        return False


def get_gff(base_url, username, password, genes, out_dir):
    features = list()

    print(f"Get gff for {len(genes)} genes")
    for key in genes:
        features.append({"uniquename": key})

    webservice_data = {"username": username, "password": password, "features": features}
    url = urllib.parse.urljoin(base_url, "annotationEditor/getGff3")

    response = requests.post(url, json=webservice_data)
    time_stamp = str(datetime.datetime.now().date())
    file_name = out_dir + "apollo_" + time_stamp + ".gff"
    if response.status_code == requests.codes.ok:
        file_handle = open(file_name, "w")
        file_handle.write(response.text)
        return file_name
    else:
        return False


def download_gff(base_url, username, password, organism, out_dir):
    webservice_data = {
        "username": username,
        "password": password,
        "type": "GFF3",
        "seqType": "genomic",
        "organism": organism,
        "output": "text",
        "exportAllSequences": "true",
        "exportGff3Fasta": "false",
    }

    url = urllib.parse.urljoin(base_url, "/IOService/write")
    response = requests.post(url, json=webservice_data)

    time_stamp = str(datetime.datetime.now().date())
    file_name = out_dir + "apollo_" + time_stamp + ".gff"

    if response.status_code == requests.codes.ok:
        file_handle = open(file_name, "w")
        file_handle.write(response.text)
        return file_name
    else:
        print(organism + ": " + str(response.text) + "\n")
        return False


def validate_gff(base_url, username, password, gff_file_path, gene_organism, moderator):
    """Validate a gff

    Returns None if there are no errors.
    """
    gff_file_object = gff_file.HandleGFF(gff_file_path, gene_organism, moderator)
    gff_file_object.read_gff_file()
    gff_file_object.scan_gff_for_errors()
    gff_file_object.scan_mrna_sequence(
        base_url=base_url, username=username, password=password
    )

    if gff_file_object.errors != {}:
        return gff_file_object
    else:
        return None


def _group_errors(errors: List[ValidationError], attrib_name: str) -> Dict[str, List[ValidationError]]:
    groups = dict()

    for error in errors:
        attrib = error.__getattribute__(attrib_name)
        if attrib in groups:
            groups[attrib].append(error)
        else:
            groups[attrib] = [error]

    return groups


def write_email_texts(errors_dict: Dict[str, List[ValidationError]], out_dir: Path) -> None:
    # First group by owner
    errors = list(errors_dict.values())
    ownership = _group_errors(errors, "owner")

    # Then write one file for each owner
    for owner, owner_errors in ownership.items():
        owner_text = _get_owner_error_text(owner, owner_errors)
        time_stamp = str(datetime.datetime.now().date())
        owner_file = out_dir / f"{owner}_{time_stamp}.error"
        with owner_file.open('w') as owner_fh:
            owner_fh.write(owner_text + "\n")


def _get_owner_error_text(owner: str, errors: List[ValidationError]) -> str:

    lines = [
        owner,
        f"Dear Annotator ({owner}),",
        (
            "***  If you've done functional annotation and not structural annotation,"
            " please ignore this message. ***"
        ),
        (
            "If you have already amended this gene annotation and think that it is correct,"
            " please ignore this message or contact us at: help@veupathdb.org to discuss any concerns."
        ),
        (
            "There is a gene annotation attributed to your account"
            " that has been edited by yourself or another annotator that currently has errors."
            " If you made any edit to this gene in the last 24 hours"
            " could you please check that the gene is correct."
        ),
        "",
    ]

    organism_groups = _group_errors(errors, "organism_name")

    for organism, organism_errors in organism_groups.items():
        organism_text = _get_organism_error_text(organism, organism_errors)
        lines.append(organism_text)

    error_text = "\n".join(lines)
    return error_text


def _get_organism_error_text(organism_name: str, errors: List[ValidationError]) -> str:

    lines = [
        f"Species: {organism_name}"
    ]

    gene_groups = _group_errors(errors, "gene_id")

    for gene_id, gene_errors in gene_groups.items():
        gene_text = _get_gene_error_text(gene_id, gene_errors)
        lines.append(gene_text)

    error_text = "\n".join(lines)
    return error_text


def _get_gene_error_text(gene_id: str, errors: List[ValidationError]) -> str:

    one_error = errors[0]
    gene_name = one_error.gene_name
    locus = one_error.locus
    lines = [
        f"Gene: {gene_name} (ID:{gene_id})\nLocation: {locus}"
    ]
    for error in errors:
        if error.mrna_id is None:
            for gene_gff_error in error.gff_error_text():
                lines.append(gene_gff_error)

    mrna_groups = _group_errors(errors, "mrna_id")

    for mrna, mrna_errors in mrna_groups.items():
        mrna_text = _get_mrna_error_text(mrna, mrna_errors)
        lines.append(mrna_text)

    error_text = "\n".join(lines)
    return error_text


def _get_mrna_error_text(mrna_id: str, errors: List[ValidationError]) -> str:

    lines = []

    for error in errors:
        for mrna_gff_error in error.gff_error_text():
            lines.append(mrna_gff_error + '.')

        for mrna_seq_error in error.sequence_error_text():
            lines.append(mrna_seq_error + '.')

    error_text = "\n".join(lines)
    return error_text


def write_summary_text(summary: AnnotatorSummary, out_dir: Path) -> None:
    """Write an email summary and a list of unfinished genes to files for an annotator.

    The files are not created if there are no annotations for that user.
    """
    owner = summary.email

    # Do not create (and so do not send) an email if there is nothing for this annotator
    if not summary.has_changes():
        print(f"No changes for {owner}")
        return

    # Create the email summary file
    file_name = out_dir / f"{owner}.summary"
    with file_name.open("w") as file_handle:
        file_handle.write(owner + "\n")
        stats = {
            "Finished protein-coding genes": summary.finished_mrna_count,
            "Unfinished protein-coding genes": summary.total_mrna_count
            - summary.finished_mrna_count,
            "Finished ncRNAs genes": summary.finished_ncrna_count,
            "Unfinished ncRNAs genes": summary.total_ncrna_count
            - summary.finished_ncrna_count,
            "Finished pseudogenes": summary.finished_pseudogene_count,
            "Unfinished pseudogenes": summary.total_pseudogene_count
            - summary.finished_pseudogene_count,
            "Non Canonical splice site": summary.non_canonical_count,
        }
        file_handle.write("Dear Annotator (" + owner + ")," + "\n")
        file_handle.write(
            "Here is a summary of your annotation in Apollo hosted at VEuPathDB.org."
            + "\n"
        )
        for item, count in stats.items():
            if count > 0:
                file_handle.write(f"{item}: {count}\n")

    # Create the gene_list file
    gene_list_name = out_dir / f"{owner}.gene_list"
    gene_list_name.touch()
    if summary.has_unfinished():
        gene_list = summary.get_unfinished()
        with open(gene_list_name, "w") as gene_list_handle:
            gene_list_handle.write("\n".join(gene_list) + "\n")


def send_email_mailgun(
    url,
    api_key,
    from_address,
    email_address,
    moderator_email_address,
    subject,
    message,
    file_attached=None,
):
    api_auth = ("api", api_key)
    mail_data = {
        "from": from_address,
        "to": email_address,
        "bcc": moderator_email_address,
        "subject": subject,
        "text": message,
    }
    req = None
    if _can_attach_file(file_attached):
        try:
            files_data = [
                (
                    "attachment",
                    ("unfinished_genes.txt", open(file_attached, "rb").read()),
                )
            ]
            req = requests.post(
                url,
                auth=api_auth,
                data=mail_data,
                files=files_data,
            )
        except ConnectionError as error:
            print(f"Connection error sending to {email_address}: {error}")
    else:
        try:
            req = requests.post(
                url,
                auth=api_auth,
                data=mail_data,
            )
        except ConnectionError as error:
            print(
                f"Connection error sending to {email_address} with file {file_attached}: {error}"
            )

    return req


def _can_attach_file(file_attached: str) -> bool:
    return (
        file_attached
        and Path(file_attached).exists()
        and Path(file_attached).stat().st_size > 0
    )
