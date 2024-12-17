
import os
import pandas as pd
import wget
import argparse
import synapseclient


def main():
    ##current AUC values are here: https://aacr.figshare.com/ndownloader/files/39996295 tabs 2 and 3
    parser = argparse.ArgumentParser()

    rawdata = 'syn64333325'
