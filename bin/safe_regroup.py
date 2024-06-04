#!/usr/bin/env python

import sys
import subprocess
import numpy as np
from biom import Table, load_table
from biom.util import biom_open

def execute_humann_regroup_table(gf_biom, group, output_file):
    """
    Written by chatGPT, thanks
    """
    # Construct the output file name
    # Construct the command as a list of arguments
    command = [
        "humann_regroup_table",
        "-i", gf_biom,
        "-g", group,
        "-o", output_file
    ]
    
    # Execute the command
    result = subprocess.run(command, capture_output=True, text=True)
    
    # Check if the command was successful
    if result.returncode == 0:
        print("Command executed successfully.")
        print(result.stdout)
    else:
        print("Command failed.")
        print(result.stderr)
    
    return result.returncode, result.stdout, result.stderr


def split_biom(original, max_samples=100):
    samples = original.ids(axis='sample')

    f = lambda id_, _: int(np.floor(list(samples).index(id_) / max_samples))

    splits = original.partition(f)

    return(splits)


def join_biom_files(input_files):

    split_bioms = []

    for x in input_files:
        split_bioms.append(load_table(x))

    base = split_bioms.pop(0)

    joined_biom = base.concat(split_bioms)

    return(joined_biom)

def main():
    args = sys.argv

    biom_fp = args[1]
    group = args[2]
    output_fp = args[3]
    if len(args) == 5:
        max_s = args[3]
    elif len(args) == 4:
        max_s = 100
    else:
        raise ValueError('Must have three or four arguments')

    print('Loading input file')
    biom_orig = load_table(biom_fp)

    print('Partitioning input file')
    biom_splits = split_biom(biom_orig, max_s)

    split_fps = []

    i = 0
    for b, t in biom_splits:
        i = i + 1
        temp_name = 'split_%s.biom' % i
        proc_name = 'split_{0}_{1}.biom'.format(i, group)
        print('Regrouping split %s' % i)
        with open(temp_name, 'w') as f:
            f.write(t.to_json('splits'))
        split_fps.append(proc_name)
        execute_humann_regroup_table(temp_name,
                                     group,
                                     proc_name)

    print('Joining split processed tables')
    joined = join_biom_files(split_fps)

    print('Saving joined table')
    with biom_open(output_fp, 'w') as f:
        tab.to_hdf5(f, 'CuratedMetagenomicData') 



if __name__ == "__main__":
    main()