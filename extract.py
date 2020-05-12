from Bio import PDB
import os
import numpy as np
import json
import functools
import argparse


cif = "3j9m"


def load_cif_from_file():

    parser = argparse.ArgumentParser()
    parser.add_argument('id', type=str,)
    parser.add_argument("--opt", type=int)
    args = parser.parse_args()
    pdbid = str.upper(args.id)

    directory = "./../cif_models/"
    filepath = directory + pdbid + ".cif"

    cifparser = PDB.FastMMCIFParser(QUIET=True)

    try:
        with open(filepath) as infile:
            structure = cifparser.get_structure(pdbid, infile)
    except:
        print("Failed to open {}".format(filepath))

    structure_atoms = PDB.Selection.unfold_entities(structure, "A")
    structure_coords = [atom.get_coord() for atom in structure_atoms]

    allchains = [chain for chain in structure.get_chains()]
    chainids = [chain.get_id() for chain in allchains]

    # chain names and  atom objects contained therein
    chain_atoms = [
        (chain.get_id(), [* chain.get_atoms()])
        for chain in allchains
    ]

    chain_coords = {
        kvpair[0]: kvpair[1] for kvpair in
        list(map(lambda chainAtomsTuple: (chainAtomsTuple[0], list(
            map(lambda atom: atom.get_coord(), chainAtomsTuple[1]))), chain_atoms))
    }

    for chain in chainids:
        if chain not in [* chain_coords.keys()]:
            print("{} is missing!".format(chain))


load_cif_from_file()
