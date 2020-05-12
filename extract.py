from Bio import PDB
import os
import numpy as np
import json
import argparse


cif = "3j9m"


def load_cif_from_file():
    parser = argparse.ArgumentParser()
    parser.add_argument('id', type=str,)
    parser.add_argument("--opt", type=int)
    args = parser.parse_args()
    print(args)


load_cif_from_file()

