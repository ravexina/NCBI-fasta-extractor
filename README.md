# NCBI fasta extractor

An script to extracts and save FASTA files from NCBI based on a query. Also generating a csv file containg details about the gathered data.
It has been designed to extract and query `nucleotide` databaes.


### How to use it

Install `biopython` using pip. Run `./main.py "query"`.

#### Example:

```text
$ ./main.py "foo"
Searching for foo

[ 7260 ] Result(s) has been found.
None of these results are in our database. All are new!

Should I fetch them? [y/N] 
```

Press [y] and it starts fetching the data. It creates a folder named `fasta` and a file named `dataset.csv`.

Dataset contain these fields:

```
key,strain,organism,isolation_source,country,year
```



