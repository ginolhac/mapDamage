#!/usr/bin/env python

import string

# from Martin Kircher, to complement DNA
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv')

letters = ('A','C','G','T','Tot')
mutations = ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G',\
        'T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') 
header=letters+mutations


