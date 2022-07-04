## IMPORTS ###

import os
import random
import csv
from generic import run_in_parallel
import numpy as np

### SOURCE ###

source = './DNMs/decode_DNMs.tsv'
chromsomes_source = './Chromosomes/WGS/Unzipped/Homo_sapiens.GRCh38.dna.chromosome.X.fa'
out_source = './CSV/trinucleotide_matrix.tsv'

### MAIN ###

def main():

    #Get TSV lines
    raw_file = open(os.path.abspath(source)).read()
    rawSplit = raw_file.split('\n')[1:]

    #First, filter to remove indels
    mutations = []
    for i in rawSplit[1:]:
        if i != '':
            line = i.split('\t')
            ref = line[2]
            alt = line[3]
            if len(ref) == 1 and len(alt) == 1:
                mutations.append(i.replace('chr', '').split('\t'))

    print (len(mutations), 'mutations')

    #Now iterate through mutations and infer trinucleotide changes
    dicts = []

    workers = int(os.cpu_count()) - 1
    processes = run_in_parallel(mutations, ['foo'], function, workers=workers)

    total_success = 0
    total_fail = 0

    for process in processes:
        output, success, fail = process.get()
        dicts.append(output)
        total_success += success
        total_fail += fail

    total_changes_dict = {}
    for d in dicts:
        for k, v in d.items():
            if k in total_changes_dict:
                total_changes_dict[k] += v
            else:
                total_changes_dict[k] = v

    #Convert dictionary to nested list
    nested_results = [[k.split('_')[0].lower(), k.split('_')[1].lower(), str(v)] for k,v in total_changes_dict.items()]
    test = np.sum([int(i[2]) for i in nested_results])
    print (test)
    print (len(nested_results))
    print (total_success, total_fail)

    #Output file
    filepath = os.path.abspath(out_source)
    headers = ['ref', 'alt', 'count']

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(i for i in headers)
        for j in nested_results:
            writer.writerow(j)

def function(mutations):

    stops = ['TAA', 'TGA', 'TAG']

    changes = {}
    success = 0
    fail = 0

    for mutation in mutations:

        chr = mutation[0]
        location = int(mutation[1])
        ref = mutation[2]
        alt = mutation[3]

        #load the relevant chromosome
        chromosome_file = open(os.path.abspath(chromsomes_source.replace('X', chr))).read()
        sequence = "".join(chromosome_file.split('\n')[1:])
        nucleotide_check = sequence[location-1]
        if ref == nucleotide_check:
            success+= 1
            trinucleotide_1_ref = sequence[location-3] + sequence[location-2] + ref
            trinucleotide_1_alt = sequence[location-3] + sequence[location-2] + alt
            change1 = trinucleotide_1_ref + '_' + trinucleotide_1_alt
            trinucleotide_2_ref = sequence[location-2] + ref + sequence[location]
            trinucleotide_2_alt = sequence[location-2] + alt + sequence[location]
            change2 = trinucleotide_2_ref + '_' + trinucleotide_2_alt
            trinucleotide_3_ref = ref + sequence[location] + sequence[location+1]
            trinucleotide_3_alt = alt + sequence[location] + sequence[location+1]
            change3 = trinucleotide_3_ref + '_' + trinucleotide_3_alt

            #Record change 1
            if trinucleotide_1_ref in stops and trinucleotide_1_alt in stops:
                if change1 not in changes:
                    changes[change1] = 1
                else:
                    changes[change1] += 1

            #Record change 2
            if trinucleotide_2_ref in stops and trinucleotide_2_alt in stops:
                if change2 not in changes:
                    changes[change2] = 1
                else:
                    changes[change2] += 1

            #Record change 3
            if trinucleotide_3_ref in stops and trinucleotide_3_alt in stops:
                if change3 not in changes:
                    changes[change3] = 1
                else:
                    changes[change3] += 1

            print (trinucleotide_1_ref, 'to', trinucleotide_1_alt, '+', trinucleotide_2_ref, 'to', trinucleotide_2_alt, '+', trinucleotide_3_ref, 'to', trinucleotide_3_alt)
        else:
            fail += 1

    return changes, success, fail

### RUN ###

if __name__ == '__main__':
    main()
