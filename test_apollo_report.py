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
import os
import shutil
import configparser
import sys
from module import annotation_quality_report as report


sys.setrecursionlimit(2500)


def make_new_directory(new_dir):
    if os.path.exists(new_dir):
        shutil.rmtree(new_dir)
        os.mkdir(new_dir)
    else:
        os.mkdir(new_dir)


def prepare_recent_gff(config, days=0):
    base_url = config['APOLLO']['base_url']
    username = config['APOLLO']['username']
    password = config['APOLLO']['password']
    config_days = int(config['SETUP']['days'])
    if days == 0:
        days = config_days
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


def main():
    parser = argparse.ArgumentParser(description='Check Apollo edits for an organism')
    
    parser.add_argument('--config', type=str, required=True,
                        help='Config file')
    
    parser.add_argument('--days', type=int,
                        help='Get genes from the last X days')
    
    parser.add_argument('--password', type=str, required=True,
                        help='Apollo password')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config)
    config['APOLLO']['password'] = args.password
    recent_genes, recent_gff = prepare_recent_gff(config, args.days)
    
    gff_file_object = report.validate_gff(
        config['APOLLO']['base_url'],
        config['APOLLO']['username'],
        config['APOLLO']['password'],
        recent_gff,
        recent_genes,
        config['EMAIL']['moderator']
    )
    print(gff_file_object)


if __name__ == '__main__':
    main()
