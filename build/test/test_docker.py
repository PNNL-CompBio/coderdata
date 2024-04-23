'''
This is a test script to validate whether or not a docker image can do the required task

'''


import os
import argparse
import subprocess


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--docker',dest='docker',default='',help='Name of docker file to test')
    parser.add_argument('--samples',dest='samples',default=False,action='store_true')
    parser.add_argument('--omics',dest='omics',default=False,action='store_true')
    parser.add_argument('--drugs',dest='drugs',default=False,action='store_true')
    parser.add_argument('--exp',dest='exp',default=False,action='store_true')

    args = parser.parse_args()

    ##make sure it can work with these test files
    test_genes='test_genes.csv'
    test_samples='test_samples.csv'
    test_drugs='test_drugs.tsv'


    ##get we generate samples?


    ##can we generate omics with genes and new samples

    ##can we generate drugs with existing drugs


    ##can we generate exp
