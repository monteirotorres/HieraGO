#!/usr/bin/env python3
# Imports
###############################################################################
import os
import re
import json
import pickle

# LICENSE
###############################################################################
'''
Crescendo Pipeline: A tool for generation of homo oligomers from pdb structures
Copyright (C) 2018 Torres, P.H.M.; Malhotra, S.; Blundell, T.L.
[The University of Cambridge]

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact info:
Department Of Biochemistry
University of Cambridge
80 Tennis Court Road
Cambridge CB2 1GA
E-mail address: monteirotorres@gmail.com
'''


# Functions
###############################################################################
def read_tsv(file):
    '''
    Receives a tsv file containing the gene names in the first column and the EC
    numbers in the second column.
    In the cases where one gene has more than one EC number, EC numbers must
    be separated by semicolons.
    In the cases where one EC number has more than gene, genes must
    be separated by spaces.
    '''
    ec_dict = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.replace('\n', '')
            locus_entry = line.split('\t')[0]
            locus_list = locus_entry.split(' ')
            ec_entry = line.split('\t')[1]
            ec_list = ec_entry.split(';')
            for locus in locus_list:
                if ec_list[0] != '':
                    ec_dict[locus] = ec_list
                else:
                    ec_dict[locus] = ['Undefined']
    return ec_dict


def count_levels(ec_dict):
    '''
    Receives a dictionary where keys are gene names (perhaps ordered locus)
    and values are lists of EC numbers. If ec number is not defined, the
    vules must be a list containing a single string 'Undefined'.
    Outputs 4 dictionaries, one per level, where keys are the partial codes
    and values are counts.
    '''
    list1 = []
    list2 = []
    list3 = []
    list4 = []
    lists = [list1, list2, list3, list4]
    dict1 = {}
    dict2 = {}
    dict3 = {}
    dict4 = {}
    dicts = [dict1, dict2, dict3, dict4]
    for l, d, n in zip(lists, dicts, range(1, 5)):
        for key, entry in ec_dict.items():
            for ec in entry:
                if n > 1:
                    if ec != 'Undefined':
                        ec = '.'.join(ec.split('.')[:n])
                        l.append(ec)
                else:
                    ec = '.'.join(ec.split('.')[:n])
                    l.append(ec)
        for ec_code in set(l):
            n = 0
            for key, entry in ec_dict.items():
                for ec in entry:
                    if ec.startswith(ec_code+'.') or ec == ec_code:
                        n += 1
            d[ec_code] = n

    return dict1, dict2, dict3, dict4


def dict2json(dict1, dict2, dict3, dict4, enzyme_dict, file_out):
    json_string = '{"name":"Enzyme Comission", "children": ['
    for key1, value1 in dict1.items():
        if key1 != 'Undefined':
            if not any(key2.startswith(key1+'.') for key2, value2 in dict2.items()):
                json_string = json_string+('{"name":"'+str(enzyme_dict[key1+'.-.-.-'])+'","code":"'+str(key1+'.-.-.-')+'","query":"'+str(dict1[key1])+'", "size": 4 }')
                continue
            else:
                json_string = json_string+('{"name":"'+str(enzyme_dict[key1+'.-.-.-'])+'","code":"'+str(key1+'.-.-.-')+'","query":"'+str(dict1[key1])+'", "children": [')
            for key2, value2 in dict2.items():
                if key2.startswith(key1):
                    if not any(key3.startswith(key2+'.') for key3, value3 in dict3.items()):
                        json_string = json_string+('{"name":"'+str(enzyme_dict[key2+'.-.-'])+'","code":"'+str(key2+'.-.-')+'","query":"'+str(dict2[key2])+'", "size": 3 }')
                        continue
                    else:
                        json_string = json_string+('{"name":"'+str(enzyme_dict[key2+'.-.-'])+'","code":"'+str(key2+'.-.-')+'","query":"'+str(dict2[key2])+'", "children": [')
                    for key3, value3 in dict3.items():
                        if key3.startswith(key2):
                            if not any(key4.startswith(key3+'.') for key4, value4 in dict4.items()):
                                json_string = json_string+('{"name":"'+str(enzyme_dict[key3+'.-'])+'","code":"'+str(key3+'.-')+'","query":"'+str(dict3[key3])+'", "size": 3 }')
                                continue
                            else:
                                json_string = json_string+('{"name":"'+str(enzyme_dict[key3+'.-'])+'","code":"'+str(key3+'.-')+'","query":"'+str(dict3[key3])+'", "children": [')
                            for key4, value4 in dict4.items():
                                if key4.startswith(key3):
                                    try:
                                        json_string = json_string+('{"name":"'+str(enzyme_dict[key4])+'","code":"'+str(key4)+'","query":"'+str(dict4[key4])+'", "size": '+str(value4)+'}')
                                    except KeyError:
                                        json_string = json_string+('{"name":"'+str(key4)+'","code":"'+str(key4)+'","query":"'+str(dict4[key4])+'", "size": '+str(value4)+'}')
                            json_string = json_string+(']}')
                    json_string = json_string+(']}')
            json_string = json_string+(']}')
    json_string = json_string+(']}')
    json_string = json_string.replace('}{', '},{')
    json_string = json_string.replace(',]', ']')
    parsed = json.loads(json_string)
    with open(file_out, 'w') as f:
        f.write(json.dumps(parsed, indent=4, ensure_ascii=False))
    return parsed


def parse_enzyme_readable(enzyme, enzclass):
    enzyme_dict = {}
    for line in open(enzyme, 'r'):
        key, value = line[:2], line[5:].rstrip()
        if key == "ID":
            id = value
        elif key == "DE":
            name = value
            if name[-1] == '.':
                name = name[:-1]
            enzyme_dict[id] = name
    for line in open(enzclass, 'r'):
        if line[0].isdigit():
            line = line.replace('\n', '')
            code = line[:9].replace(' ', '')
            if code.split('.')[1] == '-':
                name1 = line[9:-1].strip()
                enzyme_dict[code] = name1
            if code.split('.')[1] != '-' and code.split('.')[2] == '-':
                name2 = line[9:-1].strip()
                enzyme_dict[code] = name1+' - '+name2
            if code.split('.')[1] != '-' and code.split('.')[2] != '-' and code.split('.')[3] == '-':
                name3 = line[9:-1].strip()
                enzyme_dict[code] = name1+' - '+name2+' - '+name3
    enzyme_dict['none'] = 'Undefined'
    return enzyme_dict


def parse_enzyme_queriable(enzyme, enzclass):
    enzyme_dict = {}
    for line in open(enzyme, 'r'):
        key, value = line[:2], line[5:].rstrip()
        if key == "ID":
            id = value
        elif key == "DE":
            name = value
            if name[-1] == '.':
                name = name[:-1].lower()
            enzyme_dict[id] = name
    for line in open(enzclass, 'r'):
        if line[0].isdigit():
            line = line.replace('\n', '')
            code = line[:9].replace(' ', '')
            if code.split('.')[1] == '-':
                name1 = line[9:-1].strip().lower()
                enzyme_dict[code] = name1
            if code.split('.')[1] != '-' and code.split('.')[2] == '-':
                name2 = line[9:-1].strip().lower()
                enzyme_dict[code] = name1+' - '+name2
            if code.split('.')[1] != '-' and code.split('.')[2] != '-' and code.split('.')[3] == '-':
                name3 = line[9:-1].strip().lower()
                enzyme_dict[code] = name1+', '+name2+', '+name3
    enzyme_dict['none'] = 'undefined'
    return enzyme_dict


def parse_enzyme_sunburst(enzyme, enzclass):
    enzyme_dict = {}
    for line in open(enzyme, 'r'):
        key, value = line[:2], line[5:].rstrip()
        if key == "ID":
            id = value
        elif key == "DE":
            name = value
            if name[-1] == '.':
                name = name[:-1]
            enzyme_dict[id] = name
    for line in open(enzclass, 'r'):
        if line[0].isdigit():
            line = line.replace('\n', '')
            code = line[:9].replace(' ', '')
            if code.split('.')[1] == '-':
                name1 = line[9:-1].strip()
                enzyme_dict[code] = name1
            if code.split('.')[1] != '-' and code.split('.')[2] == '-':
                name2 = line[9:-1].strip()
                enzyme_dict[code] = name2
            if code.split('.')[1] != '-' and code.split('.')[2] != '-' and code.split('.')[3] == '-':
                name3 = line[9:-1].strip()
                enzyme_dict[code] = name3
    enzyme_dict['none'] = 'Undefined'
    return enzyme_dict


def full_annotation(ec_dict, out_tsv):
    if os.path.isfile(out_tsv):
        os.remove(out_tsv)
    pairs_set = set()
    for gene, ec_codes in ec_dict.items():
        if ec_codes[0] == 'Undefined':
            pairs_set.add((gene, 'none'))
        else:
            for code in ec_codes:
                p = code.split('.')
                for level in range(1, 5):
                    ndashes = 4 - level
                    pdash = p[:level] + ndashes*['-']
                    pcode = '.'.join(pdash)
                    pairs_set.add((gene, pcode))
    for pair in pairs_set:
        with open(out_tsv, 'a') as f:
            f.write(pair[0]+'\t'+pair[1]+'\n')
    return pairs_set


def pickle_dict(dict, outfile):
    with open(outfile, 'wb') as p:
        pickle.dump(dict, p, protocol=pickle.HIGHEST_PROTOCOL)


# Specify files
workdir = '/home/torres/work/abscessus/db'
enzyme = os.path.join(workdir, 'enzyme.dat')                           # Enzyme ontology from ExPASy
enzclass = os.path.join(workdir, 'enzclass.txt')                       # Enzyme partial names from ExPASy
in_tsv = os.path.join(workdir, 'uniprot_ec.tsv')                       # TSV file (obtained from UniProt)
out_json = os.path.join(workdir, 'ec.json')                            # JSON output file to create sunburst
out_tsv = os.path.join(workdir, 'ec_full_annotation.tsv')              # TSV file containing full annotation
out_readable = os.path.join(workdir, 'enzyme_dict_readable.pickle')    # Output file for picled readable enzyme dictionary
out_queriable = os.path.join(workdir, 'enzyme_dict_queriable.pickle')  # Output file for picled queriable enzyme dictionary

# Create EC ontologies
enzyme_dict_sunburst = parse_enzyme_sunburst(enzyme, enzclass)
enzyme_dict_queriable = parse_enzyme_queriable(enzyme, enzclass)
enzyme_dict_readable = parse_enzyme_readable(enzyme, enzclass)

# Read the association got from uniprot
ec_dict = read_tsv(in_tsv)

# Count the number of entries per code and create sunburst data
dict1, dict2, dict3, dict4 = count_levels(ec_dict)
dict2json(dict1, dict2, dict3, dict4, enzyme_dict_sunburst, out_json)

# Generate full annotation in the form of (gene, code) - including partial codes
pairs_set = full_annotation(ec_dict, out_tsv)

# Pickle Dictionaries
pickle_dict(enzyme_dict_readable, out_readable)
pickle_dict(enzyme_dict_queriable, out_queriable)
