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
    test_genes='/tmp/test_genes.csv'
    test_samples='/tmp/test_samples.csv'
    test_drugs='/tmp/test_drugs.tsv'


    dsname=args.docker
    if dsname in ['broad_sanger_omics','broad_sanger_exp']:
        dsname='broad_sanger;'

    if args.samples:
        ##get we generate sampless
        cmd = 'docker run -v $PWD:/tmp '+args.docker+ ' sh build_samples.sh '+test_samples
        print(cmd)
        os.system(cmd)

    if args.omics:
        ##can we generate omics with genes and new samples
        cmd = 'docker run -v $PWD:/tmp '+args.docker+ ' sh build_omics.sh '+test_genes+' /tmp/'+dsname+'_samples.csv'
        print(cmd)
        os.system(cmd)

    if args.drugs:
        ##can we generate drugs with existing drugs
        cmd = 'docker run -v $PWD:/tmp '+args.docker+ ' sh build_drugs.sh '+test_drugs
        print(cmd)
        os.system(cmd)

    if args.exp:
        ##can we generate exp
        cmd = 'docker run -v $PWD:/tmp '+args.docker+ ' sh build_exp.sh /tmp/'+dsname+'_samples.csv /tmp/'+dsname+'_drugs.tsv'
        print(cmd)
        os.system(cmd)
  

if __name__=='__main__':
    main()
