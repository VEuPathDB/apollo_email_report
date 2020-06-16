import os
import configparser
from module import annotation_quality_report as report

config = configparser.ConfigParser()
config.read('/Users/mikkel/Apollo_monitoring/annotation_quality_config.txt')
base_url = config['APOLLO']['base_url']
username = config['APOLLO']['username']
password = config['APOLLO']['password']

days = int(config['SETUP']['days'])
out_dir = config['SETUP']['dir']
static_footer = config['SETUP']['static_footer']

mailgun_url = config['MAILGUN']['url']
mailgun_key = config['MAILGUN']['api_key']

from_address = config['EMAIL']['from_address']
subject = config['EMAIL']['subject']

if not os.path.exists(out_dir):
    os.mkdir(out_dir)
else:
    _, _, last_run_files = next(os.walk(out_dir), (None, None, []))
    for file in last_run_files:
        os.unlink(out_dir + file)

recent_genes = report.get_recent_genes_from_apollo(base_url, username, password, days)
gff_file_path = False
gff_file_object = False
if recent_genes:
    gff_file_path = report.get_gff(base_url, username, password, recent_genes, out_dir)
else:
    print("No genes have been changed")
if gff_file_path:
    gff_file_object = report.validate_gff(base_url, gff_file_path, recent_genes)
else:
    print("No GFF file was downloaded")

if gff_file_object:
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

    report.sort_and_write_errors(error_lookup_table, sort_order_list, 0, out_dir)


_, _, annotation_error_files = next(os.walk(out_dir), (None, None, []))

footer_fh = open(static_footer, 'r')
footer_text = footer_fh.readlines()

for file in annotation_error_files:
    if 'error' == file[-5:]:
        file_handle = open(out_dir + file, 'r')
        email_address = file_handle.readline()
        message = file_handle.readlines()
        message += footer_text
        file_handle.close()
        email_address = 'mikkel@ebi.ac.uk'  # overwrite for testing
        report.send_email_mailgun(mailgun_url, mailgun_key, from_address, email_address, subject, message)
