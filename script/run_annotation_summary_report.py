import os
import shutil
import configparser
import json
import sys
from module import annotation_quality_report as report, gff_file

sys.setrecursionlimit(2500)


def prepare_gff(config):
    gff_dir = config['SETUP']['dir']
    organism_file = config['SETUP']['organism_file']

    organism_file_handle = open(organism_file)
    master_gff_file_name = gff_dir + '/master.gff'
    master_gff_file_handle = open(master_gff_file_name, 'w')

    gene_to_organism = dict()
    for organism_line in organism_file_handle:
        organism = organism_line.rstrip()
        gff_file_path = down_load_gff(config, organism, gff_dir)

        if not gff_file_path:
            continue

        join_gff(organism, gff_file_path, master_gff_file_handle, gene_to_organism)

    master_gff_file_handle.close()
    organism_lookup_file_name = gff_dir + '/organism_lookup.json'
    with open(organism_lookup_file_name, 'w') as outfile:
        json.dump(gene_to_organism, outfile)

    return gene_to_organism, master_gff_file_name


def import_gff_from_file(config):
    gff_dir = config['SETUP']['dir']
    organism_file = config['SETUP']['local_organism_file']

    organism_file_handle = open(organism_file)
    with open(gff_dir + '/organism_lookup.json') as json_file:
        gene_to_organism_lookup = json.load(json_file)
    master_gff_file_name = gff_dir + '/master.gff'
    master_gff_file_handle = open(master_gff_file_name, 'a')

    for organism_line in organism_file_handle:
        organism = organism_line.rstrip()
        gff_file_path = gff_dir + '/' + organism + '/' + 'apollo.gff'
        join_gff(organism, gff_file_path, master_gff_file_handle, gene_to_organism_lookup)

    master_gff_file_handle.close()
    with open(gff_dir + '/organism_lookup.json', 'w') as outfile:
        json.dump(gene_to_organism_lookup, outfile)

    return gene_to_organism_lookup, master_gff_file_name


def import_master_gff(config):
    gff_dir = config['SETUP']['dir']

    with open(gff_dir + '/organism_lookup.json') as json_file:
        gene_to_organism_lookup = json.load(json_file)
    master_gff_file_name = gff_dir + '/master.gff'

    return gene_to_organism_lookup, master_gff_file_name


def down_load_gff(config, organism_name, out_dir):
    apollo_url = config['APOLLO']['base_url']
    apollo_user_name = config['APOLLO']['username']
    apollo_password = config['APOLLO']['password']

    organism_new_dir = out_dir + '/' + organism_name + '/'
    make_new_directory(organism_new_dir)
    gff_file_name = report.download_gff(apollo_url, apollo_user_name, apollo_password, organism_name, organism_new_dir)
    return gff_file_name


def make_new_directory(new_dir):
    if os.path.exists(new_dir):
        shutil.rmtree(new_dir)
        os.mkdir(new_dir)
    else:
        os.mkdir(new_dir)


def join_gff(organism, gff_file_name, master_gff_file_handle, gene_organism):

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


def prepare_summary_emails(config, master_gff_file_name, gene_organism, file_extension):
    email_dir = config['SETUP']['dir']
    static_footer = config['SETUP']['summary_static_footer']
    footer_fh = open(static_footer, 'r')
    footer_text = footer_fh.readlines()

    gff_file_object = gff_file.HandleGFF(master_gff_file_name, gene_organism)
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
    static_footer = config['SETUP']['error_static_footer']
    footer_fh = open(static_footer, 'r')
    footer_text = footer_fh.readlines()

    gff_file_object = report.validate_gff(apollo_url, master_gff_file_name, gene_organism_lookup)

    sort_error_and_write_email_body(gff_file_object, email_dir)

    messages = list()
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


def send_emails(config, email_type, list_of_emails, use_test_email='yes'):
    email_dir = config['SETUP']['dir']
    mailgun_url = config['MAILGUN']['url']
    mailgun_key = config['MAILGUN']['api_key']
    from_address = config['EMAIL']['from_address']

    if email_type == 'summary':
        subject = config['EMAIL']['summary_subject']
    else:
        subject = config['EMAIL']['error_subject']

    for email in list_of_emails:
        email_address, email_message = email

        if email_type == 'summary':
            file_attached = email_dir + '/' + email_address + '.gene_list'
        else:
            file_attached = None

        if use_test_email == 'yes':
            email_address = config['SETUP']['test_email']  # overwrite for testing

        report.send_email_mailgun(mailgun_url, mailgun_key, from_address,
                                  email_address, subject, email_message, file_attached)


if __name__ == '__main__':
    annotation_config = configparser.ConfigParser()
    annotation_config.read('/Users/mikkel/Apollo_monitoring/annotation_summary_config.conf')

    new_gff = annotation_config['MODE']['new_gff']
    import_gff = annotation_config['MODE']['import_gff']
    run_summary = annotation_config['MODE']['run_summary']
    run_error = annotation_config['MODE']['run_error']
    new_emails = annotation_config['MODE']['new_emails']
    send_emails_flag = annotation_config['MODE']['send_emails']
    test_mode = annotation_config['MODE']['use_test_email']
    use_master_gff = annotation_config['MODE']['use_master_gff']

    if new_gff == 'yes':
        organism_lookup, master_gff = prepare_gff(annotation_config)

    if import_gff == 'yes':
        organism_lookup, master_gff = import_gff_from_file(annotation_config)

    if use_master_gff == 'yes':
        organism_lookup, master_gff = import_master_gff(annotation_config)

    if run_summary == 'yes':
        if new_emails == 'yes':
            emails = prepare_summary_emails(annotation_config, master_gff, organism_lookup, 'summary')
            if send_emails_flag == 'yes':
                send_emails(annotation_config, 'summary', emails, test_mode)

    if run_error == 'yes':
        if new_emails == 'yes':
            emails = prepare_error_emails(annotation_config, master_gff, organism_lookup,  'error')
            if send_emails_flag == 'yes':
                send_emails(annotation_config, 'error', emails, test_mode)
