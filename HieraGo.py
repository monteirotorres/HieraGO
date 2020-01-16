#!/usr/bin/env python3
# Imports
###############################################################################
import os
import pickle
import argparse
import textwrap as tw

# Description
###############################################################################
'''
This script is intended to parse data from Gene Ontology and transform it in a
hierarchical classification. It:

1 - Parses the go-basic.obo file
2 - Split into the threes different domains
3 - Defines the generations (i.e. the depth of the go-term tree)
4 - Uses the generations to codify the terms by assigning to each term the
    codes of all its parents followed by a sequential number (the number of
    the individual in that generation).
5 - Dumps a series of pickle files to be harnessed by other functionalities.

'''


# Functions
###############################################################################
def parse_go(go_base):
    '''
    Parses the terms in the 'obo' file. Does that by reading all content and
    spliting by [Term]. Outputs a dictionary in which keys are go terms
    and values are dictionaries of attributes:
    (id, name, namespace, is_obsolete, is_a = [])
    '''
    entries = open(go_base, 'r').read().split('[Term]')[1:]
    entries_dict = {}
    for entry in entries:
        attribute_dict = {}
        attributes = entry.split('\n')
        for attribute in attributes:
            if attribute.startswith('[Typedef]'):
                break
            key = str(attribute).split(': ', 1)[0]
            if key in ['id', 'name', 'namespace', 'is_obsolete']:
                value = str(attribute).split(': ', 1)[1]
            elif key == 'is_a':
                value = [str(attribute).split(': ', 1)[1].split('!')[0].replace(' ', '')]
            else:
                continue
            if key == 'is_a' and 'is_a' in attribute_dict:
                attribute_dict[key] = attribute_dict[key]+value
            else:
                attribute_dict[key] = value
        entries_dict[attribute_dict['id']] = attribute_dict
    return entries_dict


def remove_obsoletes(entries_dict):
    '''
    Simply iterates through entries dictionary and removes obsolete ones.
    '''
    for entry, data in list(entries_dict.items()):
        if 'is_obsolete' in data:
            if data['is_obsolete'] == 'true':
                del entries_dict[entry]


def split_go(entries_dict):
    '''
    Splits the parsed dictionary into the three diferent domains and assigns
    the fist letter (B, M or C) as the first part of the code.
    '''
    bio_proc = {}
    mol_func = {}
    cel_comp = {}
    rubbish = {}
    for entry, data in entries_dict.items():
        if data['namespace'] == 'biological_process':
            if 'is_a' not in data:
                data['code'] = ['B']
            else:
                data['code'] = []
            bio_proc[data['id']] = data
        elif data['namespace'] == 'molecular_function':
            if 'is_a' not in data:
                data['code'] = ['M']
            else:
                data['code'] = []
            mol_func[data['id']] = data
        elif data['namespace'] == 'cellular_component':
            if 'is_a' not in data:
                data['code'] = ['C']
            else:
                data['code'] = []
            cel_comp[data['id']] = data
        else:
            rubbish[data['id']] = data
    return bio_proc, mol_func, cel_comp, rubbish


def get_root(dict):
    '''
    From a input dictionary, retrieves the root, i.e. the term that has no
    'is_a' entry (no parents).
    '''
    for entry, data in dict.items():
        if 'is_a' not in data:
            return [data['id']]


def get_level(dict, parents_list):
    '''
    Looks for all the children of a given set of parents and return them
    as a list.
    '''
    children = []
    for entry, data in dict.items():
        for parent in parents_list:
            if 'is_a' in data:
                if parent in data['is_a']:
                    children.append(data['id'])
    return children


def get_generations(dict):
    '''
    Creates a dictionary containing the children of each generation. I.e.
    gets list of children for a set of parents, records them in children_dcit,
    makes them parents and repeat (arbitrarily until generation 20)
    '''
    children_dict = {'children0': get_root(dict)}
    for generation in range(1, 21):
        children = 'children'+str(generation)
        lastgen = 'children'+str(generation-1)
        parents = children_dict[lastgen]
        children_dict[children] = get_level(dict, parents)
    return children_dict


def codify_terms(godict, generations):
    '''
    Cycles through all individuals of each generation and assigns a code to
    them, consisting of the code of the parent followed by a sequential number.
    Returns a dictionary where keys are the unique codes and values
    are GO terms.
    '''
    code_dict = {}
    godict_cp = godict.copy()
    for g in range(0, 20):
        gen = 'children'+str(g)
        for generation, terms in generations.items():
            if generation == gen:
                n = 0
                for term in terms:
                    n += 1
                    if 'is_a' in godict_cp[term]:
                        for parent in godict_cp[term]['is_a']:
                            for parent_code in godict_cp[parent]['code']:
                                if not any(child_codes.startswith(parent_code) for child_codes in godict_cp[term]['code']):
                                    child_code = str(parent_code+'.'+str(n))
                                    godict_cp[term]['code'].append(child_code)
                        for term_code in godict_cp[term]['code']:
                            code_dict[term_code] = godict_cp[term]['id']
    return code_dict


def pickle_dict(dict, name, outpath):
    filename = os.path.join(outpath, name+'.pickle')
    with open(filename, 'wb') as p:
        pickle.dump(dict, p, protocol=pickle.HIGHEST_PROTOCOL)


def build_gochildren(hierago):
    '''
    Gets the hierarchical dictionary of GO terms and outputs a new dictionary
    in which each GO term is a key and the value is a complete list of its
    children.

    Useful to select list of terms under the selected term to query.
    '''
    gochildren = {}
    for code1, term1 in hierago.items():
        if term1 not in gochildren:
            gochildren[term1] = set()
        for code2, term2 in hierago.items():
            if code2.startswith(code1+'.') or code2 == code1:
                gochildren[term1].add(term2)
    return gochildren


def build_goparents(hierago):
    '''
    Gets the hierarchical dictionary of GO terms and outputs a new dictionary
    in which each GO term is a key and the value is a complete list of its
    parents.

    Useful to assign to a term (and ultimately to genes) a list of everything
    that it is!
    '''
    goparents = {}
    for code1, term1 in hierago.items():
        if term1 not in goparents:
            goparents[term1] = set()
        for code2, term2 in hierago.items():
            if code1.startswith(code2+'.') or code2 == code1:
                goparents[term1].add(term2)
    return goparents


def _print_generations(gendict):
    '''
    Make data for histogram of generations, for the sake of evaluation of
    number and size of generations.
    '''
    for generation, children in gendict.items():
        print(generation+': '+str(len(children)))


def _code2go(code):
    '''
    Translates HieraGO code to GO number and prints information
    from attributes dictionary
    '''
    if code.startswith('M'):
        godict = mol_func
        hierago_dict = hierago_molfunc
    if code.startswith('B'):
        godict = bio_proc
        hierago_dict = hierago_bioproc
    if code.startswith('C'):
        godict = cel_comp
        hierago_dict = hierago_celcomp
    return godict[hierago_dict[code]]


def _go2code(go):
    for key, value in hierago_molfunc.items():
        if value == go:
            print(key+': '+value)


def _len_set(hierago):
    golist = []
    for key, value in hierago.items():
        golist.append(value)
    print(name)
    print('Number of codes: '+str(len(golist)))
    print('Number of terms: '+str(len(set(golist)))+'\n')


# Main Function
###############################################################################
def main():
    print(tw.dedent("""
     ###################################################################
     ############################ HIERAGO ##############################
     ###################################################################

             Copyright (C) 2018  Torres, P.H.M.; Blundell, T.L.
                      [The University of Cambridge]

               This program comes with ABSOLUTELY NO WARRANTY


    This script is intended to parse data from Gene Ontology and
    transform it in a hierarchical classification. It:

    1 - Parses the go-basic.obo file

    2 - Split into the threes different domains

    3 - Defines the generations (i.e. the depth of the go-term tree)

    4 - Uses the generations to codify the terms by assigning to each
        term the codes of all its parents followed by a sequential
        number (the number of the individual in that generation).

    5 - Dumps a series of pickle files to be harnessed by other
        functionalities.


    Running it requires a GO classification file from the website:
    [ http://geneontology.org/page/download-ontology ]

    Also reuires an output folder where pickle files will be written
    and the selection of a classification domain:

    1 - Molecular Function
    2 - Biological Process
    3 - Cellular Component

     ###################################################################

      """))

    parser = argparse.ArgumentParser()

    parser.add_argument(dest='go_base',
                        metavar='go_base',
                        help='File containing GO rules.\n')

    parser.add_argument('-o', '--outpath',
                        dest='outpath',
                        type=str,
                        metavar='/path/to/output/files',
                        default=None,
                        help='Full path to where pickle files should be written.\n')

    parser.add_argument('-d', '--domain',
                        dest='domain',
                        type=int, default=None,
                        metavar='int',
                        help='Domain to be processed. 1, 2 or 3.\n')

    args = parser.parse_args()
    go_base = args.go_base
    outpath = args.outpath
    domain = args.domain

    assert outpath is not None, 'You must provide path where output files should go!'
    assert domain is not None, 'You must choose a domain! Read the description.'
    go_base = '/home/torres/work/go/go-basic.obo'
    entries_dict = parse_go(go_base)
    remove_obsoletes(entries_dict)
    bio_proc, mol_func, cel_comp, rubbish = split_go(entries_dict)

    if domain == 1:
        print('\n\nCalculating generations in Molecular Function GO terms...\n')
        molfuncgen = get_generations(mol_func)
        print('\nCodifying terms...\n')
        hierago_molfunc = codify_terms(mol_func, molfuncgen)
        pickle_dict(mol_func, 'molfunc', outpath)
        pickle_dict(hierago_molfunc, 'hierago_molfunc', outpath)
        print('\nBuilding parents dictionary...\n')
        goparents_molfunc = build_goparents(hierago_molfunc)
        pickle_dict(goparents_molfunc, 'goparents_molfunc', outpath)
    elif domain == 2:
        print('\n\nCalculating generations in Biological Process GO terms...\n')
        bioprocgen = get_generations(bio_proc)
        print('\nCodifying terms...\n')
        hierago_bioproc = codify_terms(bio_proc, bioprocgen)
        pickle_dict(bio_proc, 'bioproc', outpath)
        pickle_dict(hierago_bioproc, 'hierago_bioproc', outpath)
        print('\nBuilding parents dictionary...\n')
        goparents_bioproc = build_goparents(hierago_bioproc)
        pickle_dict(goparents_bioproc, 'goparents_bioproc', outpath)
    elif domain == 3:
        print('\n\nCalculating generations in Cellular Component GO terms...\n')
        celcompgen = get_generations(cel_comp)
        print('\nCodifying terms...\n')
        hierago_celcomp = codify_terms(cel_comp, celcompgen)
        pickle_dict(cel_comp, 'celcomp', outpath)
        pickle_dict(hierago_celcomp, 'hierago_celcomp', outpath)
        print('\nBuilding parents dictionary...\n')
        goparents_celcomp = build_goparents(hierago_celcomp)
        pickle_dict(goparents_celcomp, 'goparents_celcomp', outpath)

    else:
        print('\nPlease, chose a valid domain.\n')


# Execute
###############################################################################
if __name__ == '__main__':
    main()
