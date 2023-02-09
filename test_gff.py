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
import argparse
import configparser
from pathlib import Path
import sys
from module.apollo_reporter import ApolloReporter
from module import gff_file
from module import annotation_quality_report as report

sys.setrecursionlimit(2500)

moderator_name = 'moderator_name'

def main():
    parser = argparse.ArgumentParser(description='Check a single Apollo gff')
    
    parser.add_argument('--gff', type=str, required=True,
                        help='GFF to validate')
    parser.add_argument('--fasta', type=str,
                        help='Fasta with 1 CDS sequence that goes with the GFF, assuming it has a single CDS')
    parser.add_argument('--out_dir', type=str, required=True,
                        help='Outdir')
    parser.add_argument('--type', choices=('summary', 'error'), help="What do check")
    
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.type == 'error':
        gff_file_object = gff_file.HandleGFF(args.gff, {}, moderator_name)
        gff_file_object.read_gff_file()
        gff_file_object.scan_gff_for_errors()
        if args.fasta:
            gff_file_object.scan_mrna_sequence(fasta_file=args.fasta)
        
        ApolloReporter._sort_error_and_write_email_body(gff_file_object, out_dir)
        print(gff_file_object.errors)
    
    elif args.type == 'summary':
        gff_file_object = gff_file.HandleGFF(args.gff, {}, moderator_name)
        gff_file_object.read_gff_file()
        for owner, annotator_object in gff_file_object.annotators.items():
            report.write_summary_text(annotator_object, out_dir)

if __name__ == '__main__':
    main()
