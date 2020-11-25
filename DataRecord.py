import re

class DataRecord:

    def __init__(self, results):

        #  Debug
        #  print(results)

        self.key = results['GBSeq_accession-version']

        self.organism = results['GBSeq_organism']

        # Ommit month and day - Only year
        self.year = results['GBSeq_update-date'][-4:]

        # Generate country, strain, isolation source
        for rec in results['GBSeq_feature-table'][0]['GBFeature_quals']:
            key = rec['GBQualifier_name']

            if key in ['country', 'strain', 'isolation_source']:
                value = rec['GBQualifier_value']
                setattr(self, key, value)

            # Strain can also be found under isolate name
            elif key == 'isolate':
                print('strain found in isolate field')
                value = rec['GBQualifier_value']
                setattr(self, 'strain', 'isolate: ' + value)

            # else
            # it's a field that we are not intrested in

        # Set unknown fields
        for key in ['year', 'key', 'country', 'strain', 'organism', 'isolation_source']:
            if not hasattr(self, key):
                setattr(self, key, 'unknown')

        #  Remove country from data like 'India : heyderabad'
        self.country = self.country.split(':')[0]

        # Try to extract strain from title when it's not in sources
        if self.strain == 'unknown':
            self.defination = results['GBSeq_definition']
            self.extract_strain()


    def extract_strain(self):
        # Extract strain from defination if it's not in
        # source /strain or source /isolate
        p_strain = re.compile("(?<=strain\s).+?(?=(\s[^0-9]|,|\.[^0-9]))")
        match = p_strain.search(self.defination)
        if match is not None:
            print('Resolved strain from [strain]defination')
            self.strain = 'DEF: ' + match.group()
        else:
            p_isolate = re.compile("(?<=isolate\s)([A-Za-z]*\s)?.+?(?=(\s|\.[^0-9]|,))")
            match = p_isolate.search(self.defination)
            if match is not None:
                print('Resolved strain from [isolate]defination')
                self.strain = 'DEF: isolate: ' + match.group()
        # else
        # let it be unknown
