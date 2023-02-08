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
import configparser
import sys
from module.apollo_reporter import ApolloReporter

sys.setrecursionlimit(15000)

if __name__ == '__main__':
    # allow passing config file via cli arg
    if len(sys.argv) == 1:
        config_file = './config/apollo_report_config.conf'
    else:
        config_file = sys.argv[1]

    report_config = configparser.ConfigParser()
    report_config.read(config_file)

    apollo_reporter = ApolloReporter(report_config)

    # Summary annotations
    if report_config['PIPELINE']['summary_annotation'] == 'yes':
        gene_to_organism, master_gff = apollo_reporter.prepare_summary_gff()
        
        emails = apollo_reporter.prepare_summary_emails(master_gff, gene_to_organism, 'summary')
        apollo_reporter.send_emails('summary', emails)

        error_emails = apollo_reporter.prepare_error_emails(master_gff, gene_to_organism, 'error')
        apollo_reporter.send_emails('error', error_emails)
    
    # Recent annotations
    if report_config['PIPELINE']['recent_annotation'] == 'yes':
        recent_genes, recent_gff = apollo_reporter.prepare_recent_gff()
        # Write error emails
        error_emails = apollo_reporter.prepare_error_emails(recent_gff, recent_genes, 'error')
        apollo_reporter.send_emails('error', error_emails)
