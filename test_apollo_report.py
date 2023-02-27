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
import sys
from module.apollo_reporter import ApolloReporter

sys.setrecursionlimit(2500)


def main():
    parser = argparse.ArgumentParser(
        description="Check Apollo edits for an organism, no email"
    )

    parser.add_argument("--config", type=str, required=True, help="Config file")

    parser.add_argument("--password", type=str, required=True, help="Apollo password")
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config)
    config["APOLLO"]["password"] = args.password
    apollo_reporter = ApolloReporter(config)

    # Summary annotations
    if config["PIPELINE"]["summary_annotation"] == "yes":
        gene_to_organism, master_gff = apollo_reporter.prepare_summary_gff()
        emails = apollo_reporter.prepare_summary_emails(
            master_gff, gene_to_organism, "summary"
        )
        error_emails = apollo_reporter.prepare_error_emails(
            master_gff, gene_to_organism, "error"
        )

    # Recent annotations
    if config["PIPELINE"]["recent_annotation"] == "yes":
        recent_genes, recent_gff = apollo_reporter.prepare_recent_gff()
        # Write error emails
        error_emails = apollo_reporter.prepare_error_emails(
            recent_gff, recent_genes, "error"
        )


if __name__ == "__main__":
    main()
