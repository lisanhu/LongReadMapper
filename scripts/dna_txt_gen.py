#! /usr/bin/python3

import random as rnd
import argparse as arg

parser = arg.ArgumentParser()
parser.add_argument('num_characters',
                    help='number of characters to be generated', type=int)

args = parser.parse_args()

print(''.join(rnd.choices('acgt', k=args.num_characters)))
