# created by Thomas Roehl on 12/20/2021

import sys
import pandas
import urllib.parse
import urllib.request

output_id_list_file_name = sys.argv[1] # name of the UniProt-NCBI pair list output file
output_fasta_file_name = sys.argv[2] # name of the output file
input_id_list_file_name = sys.argv[3] # .csv including column of UniProt_id

## read data
print('reading data file')
data = pandas.read_csv(input_id_list_file_name)

## create query list, remove duplicates
# uniprot name list must have ids separated by spaces, example:'P40925 P40926 O43175 Q9UM73 P97793'
uniprot_id_list = ''
for index in range(len(data.UniProt_id)):
    if data.UniProt_id[index] != 'NA' and uniprot_id_list.find(data.UniProt_id[index]) > 0:
        uniprot_id_list = uniprot_id_list + " " + data.UniProt_id[index]


## code from UniProt to query online database
# instructions from UniProt: keep lists < 20,000 and remove repeats
# see https://www.uniprot.org/help/api_idmapping for list of available database conversions

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ID',
'to': 'EMBL',
'format': 'tab',
'query': uniprot_id_list
}

print('Querying UniProt Retrieve/ID mapping')
data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
ncbiresult = response.decode('utf-8')


## query NCBI EFetch
# see https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_
print('Querying NCBI EFetch')

output_list_open = open(output_id_list_file_name, 'a')

output_fasta_open = open(output_fasta_file_name, 'a')

ncbilist = ncbiresult.split('\t')
uniprotlist = uniprot_id_list.split(' ')

for index in range(len(uniprotlist)):
    ncbi = ncbilist[index]
    prot = uniprotlist[index]
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=' + ncbi + '&rettype=fasta&retmode=text'
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
       response = f.read()
    fasta = response.decode('utf-8')
    mapping = prot + ',' + ncbi + '\n'
    ouptut_list_open.write(mapping)
    output_fasta_open.write(fasta) ################## may need to add a '\n' to the end of the response
