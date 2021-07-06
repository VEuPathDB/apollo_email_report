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
import os
import shutil
import configparser
import json
import sys
from module import annotation_quality_report as report, gff_file

sys.setrecursionlimit(2500)


def prepare_summary_gff(config):
    gff_dir = config['SETUP']['dir']
    organism_file = config['SETUP']['organism_file']
    make_new_directory(gff_dir)
    organism_file_handle = open(organism_file)
    master_gff_file_name = gff_dir + '/master.gff'
    master_gff_file_handle = open(master_gff_file_name, 'w')

    gene_to_organism_lookup = dict()
    for organism_line in organism_file_handle:
        organism = organism_line.rstrip()
        gff_file_path = download_organism_gff(config, organism, gff_dir)

        if not gff_file_path:
            continue

        join_gff_to_master(organism, gff_file_path, master_gff_file_handle, gene_to_organism_lookup)

    master_gff_file_handle.close()
    organism_lookup_file_name = gff_dir + '/organism_lookup.json'
    with open(organism_lookup_file_name, 'w') as outfile:
        json.dump(gene_to_organism_lookup, outfile)

    return gene_to_organism_lookup, master_gff_file_name


def download_organism_gff(config, organism_name, out_dir):
    apollo_url = config['APOLLO']['base_url']
    apollo_user_name = config['APOLLO']['username']
    apollo_password = config['APOLLO']['password']
    clean_organism_name = organism_name.replace('/', '')
    organism_new_dir = out_dir + '/' + clean_organism_name + '/'
    make_new_directory(organism_new_dir)
    gff_file_name = report.download_gff(apollo_url, apollo_user_name, apollo_password, organism_name, organism_new_dir)
    return gff_file_name


def join_gff_to_master(organism, gff_file_name, master_gff_file_handle, gene_organism):

    gff_file_handle = open(gff_file_name, 'r')
    for line in gff_file_handle:
        fields = line.rstrip().split("\t")
        if len(fields) != 9:
            continue  # skip line as not GFF
        feature_type, owner, scaffold, strand, feature_id, parent_id, name, locus, status \
            = gff_file.extract_fields_from_gff(fields)

        if feature_type == 'gene':
            gene_organism[feature_id] = organism

        master_gff_file_handle.write(line)


def make_new_directory(new_dir):
    if os.path.exists(new_dir):
        shutil.rmtree(new_dir)
        os.mkdir(new_dir)
    else:
        os.mkdir(new_dir)


def prepare_recent_gff(config):
    base_url = config['APOLLO']['base_url']
    username = config['APOLLO']['username']
    password = config['APOLLO']['password']
    days = int(config['SETUP']['days'])
    out_dir = config['SETUP']['dir']

    make_new_directory(out_dir)

    recent_apollo_genes = report.get_recent_genes_from_apollo(base_url, username, password, days)
    gff_file_path = False
    if recent_apollo_genes:
        gff_file_path = report.get_gff(base_url, username, password, recent_apollo_genes, out_dir)
    else:
        print("No genes have been changed")
        exit()
    return recent_apollo_genes, gff_file_path


def prepare_summary_emails(config, master_gff_file_name, gene_organism, file_extension):
    email_dir = config['SETUP']['dir']
    static_footer = config['EMAIL']['summary_static_footer']
    footer_fh = open(static_footer, 'r')
    footer_text = footer_fh.readlines()

    gff_file_object = gff_file.HandleGFF(master_gff_file_name, gene_organism, config['EMAIL']['moderator'])
    gff_file_object.read_gff_file()
    write_email_body(gff_file_object, email_dir)

    messages = list()
    annotation_summary_emails = collect_files(email_dir, file_extension)
    for email in annotation_summary_emails:
        email_address, email_message = compose_message(email_dir, email, footer_text)
        messages.append((email_address, email_message))
    return messages


def write_email_body(file_object, out_dir):

    for owner, annotator_object in file_object.annotators.items():
        report.write_summary_text(annotator_object, out_dir)


def collect_files(out_dir, filter_term):
    file_filtered_list = list()
    _, _, file_list = next(os.walk(out_dir), (None, None, []))
    for file_name in file_list:
        _, file_extension = os.path.splitext(file_name)
        if '.' + filter_term == file_extension:
            file_filtered_list.append(file_name)
    return file_filtered_list


def compose_message(out_dir, email_body, footer):
    email_fh = open(out_dir + '/' + email_body, 'r')
    address = email_fh.readline().rstrip()
    message = email_fh.readlines()
    message += footer
    email_fh.close()

    return address, message


def prepare_error_emails(config, master_gff_file_name, gene_organism_lookup, file_extension):
    email_dir = config['SETUP']['dir']
    apollo_url = config['APOLLO']['base_url']
    static_footer = config['EMAIL']['error_static_footer']
    footer_fh = open(static_footer, 'r')
    footer_text = footer_fh.readlines()
    messages = list()

    gff_file_object = report.validate_gff(apollo_url, config['APOLLO']['username'], config['APOLLO']['password'],
                                          master_gff_file_name, gene_organism_lookup,
                                          config['EMAIL']['moderator'])

    if not gff_file_object:
        print('No error was found')
        return messages

    sort_error_and_write_email_body(gff_file_object, email_dir)

    annotation_error_emails = collect_files(email_dir, file_extension)
    for email in annotation_error_emails:
        email_address, email_message = compose_message(email_dir, email, footer_text)
        messages.append((email_address, email_message))
    return messages


def sort_error_and_write_email_body(gff_file_object, email_dir):

    sort_order_list = ('all', 'owner', 'organism_name', 'gene_id', 'mrna_id')
    error_object_list = list()
    error_lookup_table = dict()

    for key, value in gff_file_object.errors.items():
        error_object_list.append(value)

    error_lookup_table['all'] = error_object_list
    error_lookup_table['owner'] = list()
    error_lookup_table['organism_name'] = list()
    error_lookup_table['gene_id'] = list()
    error_lookup_table['mrna_id'] = list()

    report.sort_and_write_errors(error_lookup_table, sort_order_list, 0, email_dir)


def send_emails(config, email_type, list_of_emails):
    email_dir = config['SETUP']['dir']
    mailgun_url = config['MAILGUN']['url']
    mailgun_key = config['MAILGUN']['api_key']
    from_address = config['EMAIL']['from_address']
    email_url = config['EMAIL']['base_url']
    client_id = config['EMAIL']['client_id']
    client_secret = config['EMAIL']['client_secret']

    mode = config['SETUP']['mode']

    if email_type == 'summary':
        subject = config['EMAIL']['summary_subject']
    else:
        subject = config['EMAIL']['error_subject']

    for email in list_of_emails:
        user_id, email_message = email
        email_address = report.get_email(email_url, client_id, client_secret, user_id)
        moderators = config['EMAIL']['moderator']
        email_address = email_address + ',' + moderators
        if not email_address:
            email_address = config['EMAIL']['moderator']

        if email_type == 'summary':
            file_attached = email_dir + '/' + user_id + '.gene_list'
        else:
            file_attached = None

        if mode != 'live':
            email_address = config['EMAIL']['moderator']

        report.send_email_mailgun(mailgun_url, mailgun_key, from_address,
                                  email_address, moderator_email_address, subject, email_message, file_attached)


if __name__ == '__main__':
    # allow passing config file via cli arg
    if len(sys.argv) == 1:
        config_file = './config/apollo_report_config.conf'
    else:
        config_file = sys.argv[1]

    report_config = configparser.ConfigParser()
    report_config.read(config_file)

    if report_config['PIPELINE']['summary_annotation'] == 'yes':
        gene_to_organism, master_gff = prepare_summary_gff(report_config)
        emails = prepare_summary_emails(report_config, master_gff, gene_to_organism, 'summary')
        send_emails(report_config, 'summary', emails)
        error_emails = prepare_error_emails(report_config, master_gff, gene_to_organism, 'error')
        send_emails(report_config, 'error', error_emails)
    if report_config['PIPELINE']['recent_annotation'] == 'yes':
        recent_genes, recent_gff = prepare_recent_gff(report_config)
        error_emails = prepare_error_emails(report_config, recent_gff, recent_genes, 'error')
        send_emails(report_config, 'error', error_emails)
