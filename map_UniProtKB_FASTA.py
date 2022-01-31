# created by Thomas Roehl on 12/20/2021

import sys
import time
import pandas
import urllib.parse
import urllib.request

output_id_list_file_name = sys.argv[1] # name of the UniProt-NCBI pair list output file
output_fasta_file_name = sys.argv[2] # name of the output file
input_id_list_file_name = sys.argv[3] # .csv including column of UniProt_id

## read data
print('reading data file')
datatable = pandas.read_csv(input_id_list_file_name, keep_default_na=False)

## create query list, remove duplicates
# uniprot name list must have ids separated by spaces, example:'P40925 P40926 O43175 Q9UM73 P97793'
counter = 0
first = True
uniprotidlists = ['']
uniprotidlist = ''
uniprotidlistfull = ''
for index in range(len(datatable.UniProt_id)):
    if counter >= 5000: # split list into lists of 5000 entries to avoid UniProt server errors
        counter = 0
        if first:
            uniprotidlists[0] = uniprotidlist
            first = False
        else:
            uniprotidlists.append(uniprotidlist)
        uniprotidlist = ''
    if datatable.UniProt_id[index] != 'NA' and uniprotidlistfull.find(data.UniProt_id[index]) < 0:
        if uniprotidlist == '':
            uniprotidlist = datatable.UniProt_id[index]
        else:
            uniprotidlist = uniprotidlist + " " + datatable.UniProt_id[index]
        if uniprotidlistfull = '':
            uniprotidlistfull = datatable.UniProt_id[index]
        else:
            uniprotidlistfull = uniprotidlistfull + " " + datatable.UniProt_id[index]
if first:
    uniprotidlists[0] = uniprotidlist
else:
    uniprotidlists.append(uniprotidlist)

ncbistring = ''
for sublist in uniprotidlists:
    ## code from UniProt to query online database
    # instructions from UniProt: keep lists < 20,000 and remove repeats
    # see https://www.uniprot.org/help/api_idmapping for list of available database conversions

    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': 'ID',
    'to': 'EMBL',
    'format': 'tab',
    'query': sublist
    }

    print('Querying UniProt Retrieve/ID mapping')
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
       response = f.read()
    ncbistring = response.decode('utf-8')

ncbistlist = ncbistring.split('/n')
ncbiresult = pandas.DataFrame([row.split('\t') for row in ncbistlist])


## query NCBI EFetch
# see https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_
print('Querying NCBI EFetch')

# overwrite previous data in output files
output_list_open = open(output_id_list_file_name, 'w')
output_list_open.write('UniProtKB,NCBI_Protein\n')
output_fasta_open = open(output_fasta_file_name, 'w')
output_list_open.write('')

# open output files in append mode
output_list_open = open(output_id_list_file_name, 'a')
output_fasta_open = open(output_fasta_file_name, 'a')

ncbilist = ncbiresult[1]
uniprotlist = ncbiresult[0]

counter = 0
APIkey = '' # add in an NCBI API key here if you have one
for index in range(len(uniprotlist)):
    # slow down script to avoid error messages from NCBI
    # max 3 requests/sec if no API key, max 10 requests/sec if API key included
    if APIkey == '' and counter >= 3:
        time.sleep(1)
        counter = 1
    elif counter >= 10:
        time.sleep(1)
        counter = 1
    else:
        counter = counter + 1
    
    if ncbilist[index] != 'To': # skip any header rows ('From\tTo\n') -- lists may contain several embedded in various places
        ncbi = ncbilist[index]
        prot = uniprotlist[index]
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=' + ncbi + '&rettype=fasta&retmode=text'
        if APIkey != '':
            url = url + '&api_key=' + APIkey
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as f:
           response = f.read()
        fasta = response.decode('utf-8')
        mapping = prot + ',' + ncbi + '\n'
        output_list_open.write(mapping)
        output_fasta_open.write(fasta)
