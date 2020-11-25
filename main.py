#!/usr/bin/env python3

import os
import csv
import sys
import signal
import pickle
import time
import requests
import urllib
from Bio import Entrez
from DataRecord import DataRecord

class Extractor:

    def __init__(self, term='', db='nucleotide'):

        # It's a nice thing to tell the API who we are
        Entrez.email = "Email"
        Entrez.api_key = 'API-KEY'
        Entrez.max_tries = 1
        Entrez.sleep_between_tries = 3

        # What we are searching for
        self.term = term

        # In which database
        self.db = db

        # 10.000 is the max possible results that API supports
        # But we try our Chance
        self.max_results = 100000

        # Number of results in user query
        self.results_count = 0

        # Actual ids
        self.id_list = []

        # Path to ids list
        self.ids_path = './ids'

        # Init spreedsheet
        self.dset_path = './dataset.csv'

        # A set of ids that we already have fetch their data
        self.known_ids = self._get_known_ids()

        # Crate dataset file if not exist
        self._create_dataset()

        # Create necessary dirs if not exits
        self._create_dirs()

    def _get_known_ids(self):
        # We don't have any database of known ids yet
        if not os.path.exists(self.ids_path):
            print('Initialized an empty set for ids')
            return set()
        else:
            with open(self.ids_path, 'rb') as handle:
                known_ids = pickle.load(handle)
            return known_ids

    def _save_known_ids(self):
        with open(self.ids_path, 'wb') as handle:
            pickle.dump(self.known_ids, handle)

    def _query(self):
        handle = Entrez.esearch(db=self.db, retmax=self.max_results, term=self.term)
        query = Entrez.read(handle)
        handle.close()

        # Number of results returned from query
        self.results_count = int(query["Count"])

        # First 10 uid of results
        self.id_list = list(map(int, query["IdList"]))


    def _retrieve_fasta(self, ID):

            # Generate URL
            url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview"
            url += "&db=nuccore"
            url += "&report=fasta"
            url += "&id="
            url += str(ID)
            url += "&conwithfeat=on&hide-cdd=on"

            try:
                # Download the fasta file
                return requests.get(url, allow_redirects=True, timeout=10)
            except requests.exceptions.ReadTimeout:
                print('Connection timeout while downloading fasta.')
                return False

    def _save_fasta(self, record, fasta):

            # Create a name for saving file
            name_to_save = "./fasta/" + record.key + " " + record.country + " " + record.year + ".fasta"

            open(name_to_save, 'wb').write(fasta.content)


    def _fetch(self, ID):

            # Generate a SIGALRM if the function takes more than of 10 second to donwload
            # Because timeout and max_retry does not work on Entrez.efetch
            # We are handeling the issue ourself by limiting
            # The function call time out to 10s
            signal.alarm(10)

            # We catch the above signal and raise an exception. Here we handle the raised exception by
            # Returning a False value so we Know that we are ignoring this item for now
            # Or something else went wrong with gathering this data.
            try:
                # Create an handler which returns geneBank in xml form
                handle = Entrez.efetch(db=self.db, id=ID, rettype="gb", seq_start=1, seq_stop=2, retmode="xml")

                # Parse the actual data into python objects
                results= Entrez.read(handle)[0]

                handle.close()

                # Generate a record containing results
                return DataRecord(results)
            except ConnectionError:
                # To ignore this ID for now
                return False
            except urllib.error.URLError:
                # To ignore this ID for now
                return False
            finally:
                # reset the alarm
                signal.alarm(0)

    def _create_dataset(self):
        if not os.path.exists(self.dset_path):
            print('Creating a dataset file...')
            fields = ['key', 'strain', 'organism', 'isolation_source', 'country', 'year']
            self._append_to_dataset(fields)

    def _create_dirs(self):
        dirs = ['fasta']
        for d in dirs:
            d = './' + d
            if not os.path.exists(d):
                os.makedirs(d)


    def _append_to_dataset(self, record):
        # When fields are comming in correct form from _create_dataset
        if type(record) is type([]):
            fields = record
        # When we have to create them ourself
        else:
            fields = [record.key, record.strain, record.organism, record.isolation_source, record.country, record.year]

        with open(self.dset_path, 'a') as f:
            writer = csv.writer(f)
            writer.writerow(fields)


    def _record_exists(self, ID):
        if ID in self.known_ids:
            return True
        return False


    def _confrim(self):

        # No results
        if self.results_count == 0:
            print('Your query returned no results. quitting.')
            sys.exit(1)

        #  Find number of results FROM THIS QUERY which are already in dataset
        already_in_db = set()
        if len(self.known_ids) > 0:
            already_in_db = set(self.id_list).intersection(self.known_ids)

        # There are results, should I download them?
        print()
        print("[", str(self.results_count), "] Result(s) has been found.")

        if len(already_in_db) > 0:
            print("[", len(already_in_db), "] Number of these resutls are already in database. We will ignore them.")
        else:
            print("None of these results are in our database. All are new!")

        print()
        choice = input("Should I fetch them? [y/N] ")
        print()

        # If user wants to exit without fetching
        if choice.lower() not in ['y', 'yes']:
            print('Okay... Goodbye!')
            sys.exit(0)

        print("-"*40)
        print("\tSTARTED FEATCHING!")
        print("-"*40)
        print()
        print("-"*40)

    def _job(self, ID):

        # Check if record already exists
        if self._record_exists(ID):
            print('ID:', ID, 'Already exists, Skipping...')
            print('-'*40)
            return False

        print('Fetching', ID)
        #  returns DataRecord
        record = self._fetch(ID)

        if record == False:
            print("Cant't fetch this ID for now, Skipping...)")
            print('-'*40)
            time.sleep(10)
            return False

        #  Retrieve fasta file
        print('Downloading fasta file for', ID)
        fasta = self._retrieve_fasta(ID)
        if fasta == False:
            print('Skipping:', ID)
        else:
            print('Saving...')
            self._save_fasta(record, fasta)
            print('Done!')
        print('-'*40)

        #  Append the new record to dataset
        self._append_to_dataset(record)

        #  Append current ID to the set of records we have fetched
        self.known_ids.add(ID)

    def main(self):

        # Run the query to set the "self.id_list" and "self.results_count"
        self._query()

        self._confrim()

        for ID in self.id_list:
            self._job(ID)

        self._save_known_ids()

if __name__ == "__main__":


    if len(sys.argv) > 1:
        term = sys.argv[1]
        print("Searching for", term)
    else:
        print("Define your query!")
        sys.exit(1)

    ext = Extractor(term=term)

    def SIGINT_Handler(sig, frame):
        print("\n")
        print('Quitting...')
        ext._save_known_ids()
        sys.exit(1)

    def SIGALRM_Handler(sig, frame):
        raise ConnectionError

    # Handle SIGINT
    signal.signal(signal.SIGINT, SIGINT_Handler)
    signal.signal(signal.SIGHUP, SIGINT_Handler)
    signal.signal(signal.SIGTERM, SIGINT_Handler)
    signal.signal(signal.SIGQUIT, SIGINT_Handler)

    # Handle SIGALRM
    signal.signal(signal.SIGALRM, SIGALRM_Handler)

    try:
        ext.main()
    except KeyboardInterrupt:
        SIGINT_Handler()

