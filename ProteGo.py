#!/usr/bin/env python3
# Imports
###############################################################################
import os
import re
import json
import pickle
import argparse
import textwrap as tw
from collections import OrderedDict

# Description
###############################################################################
'''
This script is intended to generate a hierarchical picture of the Gene Ontology
terms with which proteins of a proteome have been annotated.

It parses a TSV file containg two columns:
1 - Gene Names
2 - List of go terms (separated by semicolon)

Such a table can be retrieved from uniprot as TSV or even generated beforehand.

It also needs the pickle files outputted by HieraGo script.

The final output is a JSON file containing a four-level hierarchical
organization of the proteome.
'''


# Functions
###############################################################################
def read_tsv(file):
    '''
    Receives a tsv file containing the gene names in the first column and the GO
    terms in the second column.
    In the cases where one gene has more than one GO term, GO terms must
    be separated by semicolons.
    In the cases where one GO term has more than one gene, genes must
    be separated by spaces.
    '''
    go_dict = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.replace('\n', '')
            locus_entry = line.split('\t')[0]
            locus_list = locus_entry.split(' ')
            go_entry = line.split('\t')[1]
            go_list = go_entry.split('; ')
            for locus in locus_list:
                if go_list[0] != '':
                    go_dict[locus] = go_list
                else:
                    go_dict[locus] = ['Undefined']
    return go_dict


def reverse_hiearago(hierago):
    '''
    Receives a dictionary with hierarchical GO entries as keys and regular
    GO terms as values and reverses it, creating a dictionary with GO terms as
    keys and the values are lists of codes.
    '''
    hierago_reverse = {}
    for key, entry in hierago.items():
        if entry in hierago_reverse:
            hierago_reverse[entry].append(key)
        else:
            hierago_reverse[entry] = []
            hierago_reverse[entry].append(key)
    return hierago_reverse


def codify_proteome(go_dict, hierago_reverse):
    '''
    Receives the dictionary parsed from csv file and the dictionary correlating
    GO terms with hierarchical codes and returns a dictionary assigning  a list
    of hierarchical codes for each protein
    '''
    coded_godict = {}
    for gene, go_terms in go_dict.items():
        for go in go_terms:
            if go in hierago_reverse:
                if gene in coded_godict:
                    for code in hierago_reverse[go]:
                        if code not in coded_godict[gene]:
                            coded_godict[gene].append(code)
                            # coded_godict[gene].append(code+'.-')
                else:
                    coded_godict[gene] = []
                    for code in hierago_reverse[go]:
                        if code not in coded_godict[gene]:
                            coded_godict[gene].append(code)
                            # coded_godict[gene].append(code+'.-')
    for gene, go_terms in go_dict.items():
        if gene not in coded_godict:
            coded_godict[gene] = ['Undefined']
    return coded_godict


def count_levels(coded_godict):
    '''
    Receives a dictionary where keys are gene names (perhaps ordered locus)
    and values are lists of hierarchical GO codes.
    Outputs 4 dictionaries, one per level, where keys are the partial codes
    and values are counts.
    '''
    dicts = OrderedDict()
    for g in range(1, 7):
        s = set()
        d = {}
        for gene, codes in coded_godict.items():
            for code in codes:
                levels = code.split('.')[:g]
                if len(levels) == g:
                    part = '.'.join(levels)
                    s.add(part)
        for part in s:
            geneset = set()
            # n = 0
            for gene, codes in coded_godict.items():
                for code in codes:
                    if code.startswith(part+'.') or code == part:
                        # n += 1
                        geneset.add(gene)
            d[part] = len(geneset)  # n
        dicts[g] = d
    return dicts


def back_count(dicts):
    for layer, entries in reversed(dicts.items()):
        if layer > 2:
            ancestor_layer = dicts[layer-1]
            for ancestor_code, ancestor_count in ancestor_layer.items():
                n = 0
                for code, count in dicts[layer].items():
                    if code.startswith(ancestor_code+'.'):
                        n = n+count
                if ancestor_count > n:
                    remainder_count = ancestor_count - n
                    remainder_code = ancestor_code+'.-'
                    dicts[layer][remainder_code] = remainder_count
    return dicts


def dict2json(root_dict, dict1, dict2, dict3, dict4, dict5,  main_dict, hierago, json_out, domain):
    '''
    Creates the data to be read by java D3 library and plotted as a sunburst.
    Does that by iterating through the codes in the hierarchical dictionary
    checking their children by comparing the partial codes and outputting
    the data in JSON format:

    {name: AAAA, children:[
        {name: BBBB, children:[
            {name CCCC, size:X}   <----- this will be last level in sunburst
        ]}
    ]}

    '''

    if domain == 'molfunc':
        dname = 'Molecular Function'
        dcode = 'M'
        dterm = 'GO:0003674'
    elif domain == 'bioproc':
        dname = 'Biological Process'
        dcode = 'B'
        dterm = 'GO:0008150'
    elif domain == 'celcomp':
        dname = 'Cellular Component'
        dcode = 'C'
        dterm = 'GO:0005575'

    json_string = ''
    json_string = json_string+('{"name": "'+dname+'", "code": "'+dcode+'", "term": "'+dterm+'", "query":'+str(root_dict[dcode])+', "children": [')
    for key1, value1 in dict1.items():
        if not any(key2.startswith(key1+'.') for key2, value2 in dict2.items()):
            json_string = json_string+('{"name":"'+main_dict[hierago[key1]]['name']+'", "code": "'+key1+'", "term": "'+hierago[key1]+'", "query": '+str(value1)+', "size": 5 }')
            continue
        else:
            json_string = json_string+('{"name":"'+main_dict[hierago[key1]]['name']+'", "code": "'+key1+'", "term": "'+hierago[key1]+'", "query": '+str(value1)+', "children": [')
        for key2, value2 in dict2.items():
            if key2.startswith(key1+'.'):
                if not any(key3.startswith(key2+'.') for key3, value3 in dict3.items()):
                    json_string = json_string+('{"name":"'+main_dict[hierago[key2]]['name']+'", "code": "'+key2+'", "term": "'+hierago[key2]+'", "query": '+str(value2)+', "size": 4 }')
                    continue
                else:
                    json_string = json_string+('{"name":"'+main_dict[hierago[key2]]['name']+'", "code": "'+key2+'", "term": "'+hierago[key2]+'", "query": '+str(value2)+', "children": [')
                for key3, value3 in dict3.items():
                    if key3.startswith(key2+'.'):
                        if not any(key4.startswith(key3+'.') for key4, value4 in dict4.items()):
                            json_string = json_string+('{"name":"'+main_dict[hierago[key3]]['name']+'", "code": "'+key3+'", "term": "'+hierago[key3]+'", "query": '+str(value3)+', "size": 3 }')
                            continue
                        else:
                            json_string = json_string+('{"name":"'+main_dict[hierago[key3]]['name']+'", "code": "'+key3+'", "term": "'+hierago[key3]+'", "query": '+str(value3)+', "children": [')
                        for key4, value4 in dict4.items():
                            if key4.startswith(key3+'.'):
                                if not any(key5.startswith(key4+'.') for key5, value5 in dict5.items()):
                                    json_string = json_string+('{"name":"'+main_dict[hierago[key4]]['name']+'", "code": "'+key4+'", "term": "'+hierago[key4]+'", "query": '+str(value4)+', "size": 2 }')
                                    continue
                                else:
                                    json_string = json_string+('{"name":"'+main_dict[hierago[key4]]['name']+'", "code": "'+key4+'", "term": "'+hierago[key4]+'", "query": '+str(value4)+', "children": [')
                                for key5, value5 in dict5.items():
                                    if key5.startswith(key4+'.'):
                                        json_string = json_string+('{"name":"'+main_dict[hierago[key5]]['name']+'", "code": "'+key5+'", "term": "'+hierago[key5]+'", "query": '+str(value5)+', "size": 1 }')
                                json_string = json_string+(']}')
                        json_string = json_string+(']}')
                json_string = json_string+(']}')
        json_string = json_string+(']}')
    json_string = json_string+(']}')
    json_string = json_string.replace('}{', '},{')
    parsed = json.loads(json_string)
    with open(json_out, 'w') as f:
        f.write(json.dumps(parsed, indent=4))
    return parsed


def get_full_annotation(go_dict, goparents, full_out):
    gene_parents = {}
    for gene, go_list in go_dict.items():
        gofull = set()
        for go in go_list:
            if go in goparents:
                gofull.add(go)
                for parent in goparents[go]:
                    gofull.add(parent)
        gene_parents[gene] = gofull
    with open(full_out, 'w') as f:
        for gene, codes in gene_parents.items():
            for code in codes:
                f.write(gene+'\t'+code+'\n')


def _retrieve_genes(query_code):
    uniques = set()
    for gene, codes in coded_godict.items():
        for code in codes:
            if code.startswith(query_code+'.') or code == query_code:
                print(gene+': '+code)
                uniques.add(gene)
    print('Number of genes that will be queried: '+str(len(uniques)))


def _count_undefineds():
    n = 0
    for gene in coded_godict:
        if coded_godict[gene] == ['Undefined']:
            # print(gene)
            n += 1
    print(n)


def _count_defineds():
    n = 0
    for gene, codes in coded_godict.items():
        for code in codes:
            if code != 'Undefined':
                n += 1
                # break
    print(n)


def _count_startswith(query_code):
    n = 0
    for gene, codes in coded_godict.items():
        for code in codes:
            if code.startswith(query_code):
                n += 1
    print(n)


def _find_genes(term):
    full_annotation = get_full_annotation(go_dict, goparents)
    results = set()
    for gene, terms in full_annotation.items():
        if term in terms:
            results.add(gene)
    return results


# Main Function
###############################################################################
def main():
    print(tw.dedent("""
     ###################################################################
     ############################ PROTEGO ##############################
     ###################################################################

             Copyright (C) 2018  Torres, P.H.M.; Blundell, T.L.
                      [The University of Cambridge]

               This program comes with ABSOLUTELY NO WARRANTY


    This script is intended to parse GO annotation data from Uniprot and
    transform it into a hierarchical classification. It:

    1 - Parses tsv file containing gene names and GO terms.
    2 - Uses pickled HIERAGO dictionaries to codify the GO terms in
        proteome.
    3 - Defines number of indivials in each generation and their
        children.
    4 - Uses that information to build JSON file to be read by D3.
    5 - Outputs a list of genes annotated with all GO terms, not just
        terminal annotations.

    Running it requires te path to pickle files outputted from hierago.

    Also an output folder where JSON and tsv files will be written.

    Select a domain:
    1 - Molecular Function
    2 - Biological Process
    3 - Cellular Component

     ###################################################################

      """))

    parser = argparse.ArgumentParser()

    parser.add_argument(dest='input_file',
                        metavar='input_file',
                        help='CSV file containing gene names (col1) and GO terms (col2)\n')

    parser.add_argument('-p', '--picklepath',
                        dest='picklepath',
                        type=str,
                        metavar='/path/to/pickle/files',
                        default=None,
                        help='Full path from where pickle files should be read.\n')

    parser.add_argument('-o', '--outpath',
                        dest='outpath',
                        type=str,
                        metavar='/path/to/output/files',
                        default=None,
                        help='Full path to where JSON and TSV files should be written.\n')

    parser.add_argument('-d', '--domain',
                        dest='domain',
                        type=int, default=None,
                        metavar='int',
                        help='Domain to be processed. 1, 2 or 3.\n')

    args = parser.parse_args()
    input_file = args.input_file
    picklepath = args.picklepath
    outpath = args.outpath
    domain = args.domain

    assert picklepath is not None, 'You must provide path where HIERAGO pickle files are located!'
    assert outpath is not None, 'You must provide path where output files should go!'
    assert domain is not None, 'You must choose a domain! Read the description.'

    if domain == 1:
        domain = 'molfunc'
    elif domain == 2:
        domain = 'bioproc'
    elif domain == 3:
        domain = 'celcomp'
    else:
        print('\nPlease choose a valid domain!\n')

    json_out = os.path.join(outpath, domain+'.json')
    full_out = os.path.join(outpath, domain+'_full_annotation.tsv')

    main_dict_pickle = os.path.join(picklepath, domain+'.pickle')
    with open(main_dict_pickle, 'rb') as p:
        main_dict = pickle.load(p)

    hierago_pickle = os.path.join(picklepath, 'hierago_'+domain+'.pickle')
    with open(hierago_pickle, 'rb') as p:
        hierago = pickle.load(p)

    goparents_pickle = os.path.join(picklepath, 'goparents_'+domain+'.pickle')
    with open(goparents_pickle, 'rb') as p:
        goparents = pickle.load(p)

    go_dict = read_tsv(input_file)
    hierago_reverse = reverse_hiearago(hierago)
    coded_godict = codify_proteome(go_dict, hierago_reverse)
    dicts = count_levels(coded_godict)
    dict2json(dicts[1], dicts[2], dicts[3], dicts[4], dicts[5], dicts[6], main_dict, hierago, json_out, domain)
    get_full_annotation(go_dict, goparents, full_out)


# Execute
###############################################################################
if __name__ == '__main__':
    main()
