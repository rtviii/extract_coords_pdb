from Bio import PDB
import pickle
import pandas as pd
import os
import numpy as np
import json
import functools
import argparse


# PATH TO CIF FILES FROM PDB
directory = "./cif_models/"


def save_molecule_as_csv(pdbid, partition):
    for kvpair in partition:
        if not os.path.exists("{}".format(pdbid)):
            try:
                os.makedirs('{}'.format(pdbid))
            except:
                print("Failed to create directory for molecule. Exiting")

        df = pd.DataFrame(kvpair[1])
        df.columns = ['x', 'y', 'z']
        outpath = './{}/{}.csv'.format(pdbid, kvpair[0])
        df.to_csv(outpath, index=True)
        print("Saved chain {} coordinate to {}.csv ".format(kvpair[0], kvpair[0]))

# Usage: call with the id of the structure
# Ex. python3 extract_coordinates.py 5jvg
# Deposits a dictionary of subchains and their corresponding coordinates to /coordinate_dicts


def load_cif_from_file():

    # Parse the pdbid argument
    parser = argparse.ArgumentParser()
    parser.add_argument('id', type=str,)
    args = parser.parse_args()
    pdbid = str.upper(args.id)
    filepath = directory + pdbid + ".cif"

    # Initializing pdb parsermole
    cifparser = PDB.FastMMCIFParser(QUIET=True)
    try:
        with open(filepath) as infile:
            structure = cifparser.get_structure(pdbid, infile)
    except:
        print("Failed to open {}".format(filepath))

    # Decomposing the structure by chain
    structure_atoms = PDB.Selection.unfold_entities(structure, "A")
    structure_coords = [atom.get_coord() for atom in structure_atoms]
    allchains = [chain for chain in structure.get_chains()]
    chainids = [chain.get_id() for chain in allchains]

    # chain names and atom objects contained therein
    chain_atoms = [(chain.get_id(), [* chain.get_atoms()])
                   for chain in allchains]

    def atomArrToCord(array_of_atoms):
        return list(map(lambda atom:
                        atom.get_coord(), array_of_atoms)
                    )
    chainCoordinates = list( map((lambda chainAtomTuple: (chainAtomTuple[0], atomArrToCord(
        chainAtomTuple[1]))), chain_atoms) )


    chain_coordinates_object = {}

    for namecoordpair in chainCoordinates:
        chain_coordinates_object[namecoordpair[0]] =  namecoordpair[1]


    save_molecule_as_csv(pdbid, chainCoordinates)

    # # Sanity check
    # for chain in chainids:
    #     if chain not in [* chain_coords.keys()]:
    #         print("{} is missing!".format(chain))

    # if not os.path.exists('coordinate_dicts'):
    #     try:
    #         os.makedirs('coordinate_dicts')
    #     except:
    #         print("Failed to create deposition directory. Exiting")

    # picklename = "./coordinate_dicts/{}.pkl".format(pdbid)
    # output = open(picklename, 'wb')
    # pickle.dump(chain_coords, output)
    # print("Saved successfully to {}/coordinate_dicts/{}.pkl".format(os.getcwd(), pdbid))
    # output.close()


load_cif_from_file()
