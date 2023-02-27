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
import json
import sys
from module import annotation_quality_report as report

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

    organisms = report.get_organisms(
        config["APOLLO"]["base_url"],
        config["APOLLO"]["username"],
        config["APOLLO"]["password"],
    )
    print(json.dumps(organisms, indent=1))


if __name__ == "__main__":
    main()
