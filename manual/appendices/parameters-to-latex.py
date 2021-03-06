#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import collections
import re
import sys
import yaml

# Load YAML maps as ordered dictionaries to preserve the grouping of parameters
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))

yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)


def replace(s):
    replacements = [
        (r"@MATH\[([^\]]*)\]@",     r"$\1$"),
        (r"@REF\[(.*)\]@",         r"\cite{\1}"),
    ]
    result = s
    for r in replacements:
        result = re.sub(r[0], r[1], result)

    return result


def main():
    if not len(sys.argv) == 2:
        print('2 arguments required, {} given!'.format(len(sys.argv) - 1))
        exit(-1)

    infilename = sys.argv[1]

    with open(infilename) as data_file:
        try:
            parameters = yaml.load(data_file)
        except yaml.YAMLERROR as e:
            print(e)
            exit(-1)

    blacklisted_keys = ['@metadata@']
    keys = [key for key in parameters if key not in blacklisted_keys]

    description = r'NO DESCRIPTION'
    if '@metadata@' in parameters:
        if 'description' in parameters['@metadata@']:
            description = replace(parameters['@metadata@']['description'])

        print(r'{%')
        print(r'\renewcommand{\arraystretch}{1.1}')
        print(r'\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}')
        print(r'\begin{center}')
        print(r'\begin{tabular}{@{} p{.5\textwidth} p{.1\textwidth} L{.4\textwidth} @{}}')
        print(r'\toprule')
        print(r'    \textbf{Qualified Name} & \textbf{Parameter} & \textbf{Description} \\')
        print(r'\midrule')
        for key in keys:
            p = parameters[key]

            latex = '???'
            if 'latex' in p:
                latex = p['latex']

            comment = ''
            if 'comment' in p:
                comment = replace(p['comment'])

            print(r'    \verb|{0}| & {1} & {2} \\'.format(key, latex, comment))

        print(r'\bottomrule')
        print(r'\end{tabular}')
        print(r'\end{center}')
        print(r'\captionof{{table}}{{ {0} }}'.format(description))
        print(r'}')

    exit(0);

if __name__ == '__main__':
    main()
