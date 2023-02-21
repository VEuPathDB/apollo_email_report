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
import requests
import datetime
import urllib.parse
import re
from module import gff_file
from module.annotator import AnnotatorSummary


def get_recent_genes_from_apollo(base_url, username, password, days=1):

    webservice_data = {'username': username, 'password': password, 'days': days}
    url = urllib.parse.urljoin(base_url, 'annotationEditor/getRecentAnnotations')

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

    webservice_data = {'username': username, 'password': password}
    url = urllib.parse.urljoin(base_url, 'organism/findAllOrganisms')

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
    user_obj = re.match(r'.*\.(\d+?)$', user_id, flags=0)

    if user_obj:
        upenn_id = int(user_obj.group(1))
    else:
        return False

    body = {"client_id": client_id,
            "client_secret": client_secret,
            "query": {
                "userId": upenn_id
            }}
    response = requests.post(base_url, json=body)
    if response.status_code == requests.codes.ok:
        return response.json()["email"]
    else:
        return False


def get_gff(base_url, username, password, genes, out_dir):

    features = list()

    print(f"Get gff for {len(genes)} genes")
    for key in genes:
        features.append({'uniquename': key})

    webservice_data = {'username': username, 'password': password, 'features': features}
    url = urllib.parse.urljoin(base_url, 'annotationEditor/getGff3')

    response = requests.post(url, json=webservice_data)
    time_stamp = str(datetime.datetime.now().date())
    file_name = out_dir + 'apollo_' + time_stamp + '.gff'
    if response.status_code == requests.codes.ok:
        file_handle = open(file_name, 'w')
        file_handle.write(response.text)
        return file_name
    else:
        return False



def download_gff(base_url, username, password, organism, out_dir):

    webservice_data = {'username': username,
                       'password': password,
                       'type': 'GFF3',
                       'seqType': 'genomic',
                       'organism': organism,
                       'output': 'text',
                       'exportAllSequences': 'true',
                       'exportGff3Fasta': 'false'}

    url = urllib.parse.urljoin(base_url, '/IOService/write')
    response = requests.post(url, json=webservice_data)

    time_stamp = str(datetime.datetime.now().date())
    file_name = out_dir + 'apollo_' + time_stamp + '.gff'

    if response.status_code == requests.codes.ok:
        file_handle = open(file_name, 'w')
        file_handle.write(response.text)
        return file_name
    else:
        print(organism + ': ' + str(response.text) + '\n')
        return False


def validate_gff(base_url, username, password, gff_file_path, gene_organism, moderator):
    """Validate a gff
    
    Returns None if there are no errors.
    """
    gff_file_object = gff_file.HandleGFF(gff_file_path, gene_organism, moderator)
    gff_file_object.read_gff_file()
    gff_file_object.scan_gff_for_errors()
    gff_file_object.scan_mrna_sequence(base_url=base_url,  username=username, password=password)

    if gff_file_object.errors != {}:
        return gff_file_object
    else:
        return None


def sort_and_write_errors(dict_of_list, order_of_lists, index, out_dir: Path, file_handle=None):
    file_handle = file_handle
    # print("length of list:", order_of_lists[index], len(dict_of_list[order_of_lists[index]]))
    if index < 0:
        return False
    elif index == 0 and len(dict_of_list[order_of_lists[index]]) == 0:
        return False
    elif len(dict_of_list[order_of_lists[index]]) > 0 and index + 1 < len(order_of_lists):
        # print("copy to and writing out for", order_of_lists[index + 1])
        copy_function(dict_of_list, order_of_lists, index)
        file_handle = write_function(dict_of_list, order_of_lists, index + 1, out_dir, file_handle)
        sort_and_write_errors(dict_of_list, order_of_lists, index + 1, out_dir, file_handle)
    else:
        # print("clear list", order_of_lists[index])
        # print("going up to", order_of_lists[index - 1])
        dict_of_list[order_of_lists[index]].clear()
        sort_and_write_errors(dict_of_list, order_of_lists, index - 1, out_dir, file_handle)
    if file_handle:
        file_handle.close()


def copy_function(dict_of_list, order_of_lists, index):
    search_term = order_of_lists[index + 1]
    search_value = str()

    for error_object in dict_of_list[order_of_lists[index]]:
        # print("search value", search_value, "search term", error_object.__dict__[search_term])
        if not search_value:
            search_value = error_object.__dict__[search_term]
            # print("setting search value to", search_value)
            dict_of_list[order_of_lists[index + 1]].append(error_object)
            continue
        if error_object.__dict__[search_term] == search_value:
            # print("copy to", order_of_lists[index + 1])
            dict_of_list[order_of_lists[index + 1]].append(error_object)

    if dict_of_list[order_of_lists[index + 1]]:
        for delete_object in dict_of_list[order_of_lists[index + 1]]:
            # print("delete from", order_of_lists[index])
            dict_of_list[order_of_lists[index]].remove(delete_object)


def write_function(dict_of_list, order_of_list, index, out_dir: Path, file_handle=None):
    file_handle = file_handle
    if order_of_list[index] == 'owner':
        if file_handle:
            file_handle.close()
        owner = dict_of_list['owner'][0].owner
        time_stamp = str(datetime.datetime.now().date())
        file_name = out_dir / f"{owner}_{time_stamp}.error"
        file_handle = file_name.open('a')
        file_handle.write(owner + "\n")
        file_handle.write("Dear Annotator (" + owner + ")," + "\n")
        file_handle.write("***  If you've done functional annotation and not structural annotation, please ignore this message. ***" + "\n")
        file_handle.write("If you have already amended this gene annotation and think that it is correct, please ignore this message or contact us at: help@veupathdb.org to discuss any concerns." + "\n")
        file_handle.write("There is a gene annotation attributed to your account that has been edited by yourself or another annotator that currently has errors. If you made any edit to this gene in the last 24 hours could you please check that the gene is correct." + "\n")
    elif order_of_list[index] == 'organism_name':
        organism_name = dict_of_list['organism_name'][0].organism_name
        organism_str = "Species: {}\n".format(organism_name)
        file_handle.write(organism_str)
    elif order_of_list[index] == 'gene_id':
        gene_name = dict_of_list['gene_id'][0].gene_name
        gene_id = dict_of_list['gene_id'][0].gene_id
        locus = dict_of_list['gene_id'][0].locus
        gene_str = "Gene: {} (ID:{})\nLocation: {}\n".format(gene_name, gene_id, locus)
        file_handle.write(gene_str)
        for error_object in dict_of_list['gene_id']:
            if error_object.mrna_id is None:
                for string in error_object.gff_error_text():
                    file_handle.write(string)

    elif order_of_list[index] == 'mrna_id':
        for error_object in dict_of_list['mrna_id']:

            for string in error_object.gff_error_text():
                file_handle.write(string)

            for string in error_object.sequence_error_text():
                file_handle.write(string)

    return file_handle


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
    with file_name.open('w') as file_handle:
        file_handle.write(owner + "\n")
        stats = {
            "Finished protein-coding genes": summary.finished_mrna_count,
            "Unfinished protein-coding genes": summary.total_mrna_count - summary.finished_mrna_count,
            "Finished ncRNAs genes": summary.finished_ncrna_count,
            "Unfinished ncRNAs genes": summary.total_ncrna_count - summary.finished_ncrna_count,
            "Finished pseudogenes": summary.finished_pseudogene_count,
            "Unfinished pseudogenes": summary.total_pseudogene_count - summary.finished_pseudogene_count,
            "Non Canonical splice site": summary.non_canonical_count,
        }
        file_handle.write('Dear Annotator (' + owner + '),' + "\n")
        file_handle.write('Here is a summary of your annotation in Apollo hosted at VEuPathDB.org.' + "\n")
        for (item, count) in stats.items():
            if count > 0:
                file_handle.write(f"{item}: {count}\n")
    
    # Create the gene_list file
    gene_list_name = out_dir / f"{owner}.gene_list"
    gene_list_name.touch()
    if summary.has_unfinished():
        gene_list = summary.get_unfinished()
        with open(gene_list_name, 'w') as gene_list_handle:
            gene_list_handle.write("\n".join(gene_list) + "\n")


def send_email_mailgun(url, api_key, from_address, email_address, moderator_email_address, subject, message, file_attached=None):

    if file_attached:
        return requests.post(url, auth=("api", api_key), files=[("attachment", ("unfinished_genes.txt",
                                                                                open(file_attached, "rb").read()))],
                             data={"from": from_address, "to": email_address, "bcc": moderator_email_address, "subject": subject, "text": message})
    else:
        return requests.post(url, auth=("api", api_key), data={"from": from_address, "to": email_address, "bcc": moderator_email_address,
                                                               "subject": subject, "text": message})
