# apollo_email_report
A tool used to monitor and validate annotation in Apollo(genome annotation web-based editor). A report will be emailed to the annotator(owner of the annotation).
The tool has two modes recent and summary. Recent will download annotation for a given number of days and validate it, if there is any errors an email will be sent to the annotator.
Summary will download all annotation for the organisms given in 'organism_file.text' and sent a summary and any errors to the annotator. MAILGUN is used to send emails and a MAILGUN account is required.

## Installation
Setup the config file: ./config/apollo_report_config.conf

```bash
pip install -r requirements.txt
```

## Usage
```bash
python3  run_apollo_report.py
```
