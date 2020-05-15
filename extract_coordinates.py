from Bio import PDB
import pickle
import pandas as pd
import os
import numpy as np
import json
import functools
import argparse


# PATH TO CIF FILES (acquired from pdb)
directory = "./cif_models/"


def pickle_coordinates(pdbid, partition):
    picklename = "./coordinate_dicts/{}.pkl".format(pdbid)
    output     = open(picklename, 'wb')
    pickle.dump(partition, output)
    print("Saved successfully to {}/coordinate_dicts/{}.pkl".format(os.getcwd(), pdbid))
    output.close()

    
def save_molecule_as_csv(pdbid, partition):
    print("Length of partition: ", len(partition))
    print("Parsing to .csv...")


    df = pd.DataFrame({key: pd.Series(value) for key, value in partition.items()})

    for kvpair in partition:
        if not os.path.exists("{}".format(pdbid)):
            try:
                os.makedirs('{}'.format(pdbid))
            except:
                print("Failed to create directory for molecule. Exiting")
    outpath = './coordinate_csvs/{}.csv'.format(pdbid, pdbid)
    df.to_csv(outpath, index=True)
    print("Saved coordinates to {}.csv ".format(pdbid))


# Usage: call with the id of the structure
# Ex. python3 extract_coordinates.py 5jvg
# Deposits a dictionary of subchains and their corresponding coordinates to /coordinate_dicts

def load_cif_from_file():

    # Parse the pdbid argument
    parser   = argparse.ArgumentParser()
    parser.add_argument('id', type=str,)
    args     = parser.parse_args()
    pdbid    = str.upper(args.id)
    filepath = directory + pdbid + ".cif"

    # Initializing pdb parsermole
    cifparser = PDB.FastMMCIFParser(QUIET=True)
    try:
        with open(filepath) as infile:
            structure = cifparser.get_structure(pdbid, infile)
    except:
        print("Failed to open {}".format(filepath))

    # Decomposing the structure by chain
    structure_atoms  = PDB.Selection.unfold_entities(structure, "A")
    structure_coords = [atom.get_coord() for atom in structure_atoms]
    allchains        = [chain for chain in structure.get_chains()]
    chainids         = [chain.get_id() for chain in allchains]

    # chain names and atom objects contained therein
    chain_atoms = [(chain.get_id(), [* chain.get_atoms()]) for chain in allchains]

    def atomArrToCord(array_of_atoms):
        return list(map(lambda atom:
                        atom.get_coord(), array_of_atoms)
                    )
    chainCoordinates = list(map((lambda chainAtomTuple: (chainAtomTuple[0], atomArrToCord(
        chainAtomTuple[1]))), chain_atoms))

    chain_coordinates_object = {}

    for namecoordpair in chainCoordinates:
        chain_coordinates_object[namecoordpair[0]] = namecoordpair[1]
    save_molecule_as_csv(pdbid, chain_coordinates_object)


    # Sanity check
    for chain in chainids:
        if chain not in [* chain_coordinates_object.keys()]:
            print("{} is missing!".format(chain))

    if not os.path.exists('coordinate_dicts'):
        try:
            os.makedirs('coordinate_dicts')
        except:
            print("Failed to create deposition directory. Exiting")


load_cif_from_file()
