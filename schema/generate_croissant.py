'''
script that generates the croissant metadata file
'''


import json
import argparse
import yaml

def main():
    parser = argparse.ArgumentParser('Script to generate croissant metadata for CoderData')
    parser.add_argument('--figsharePath',dest='figshare',help='Path prefix to figshare to generate URLs for each file')
    parser.add_argument('--yamlFile',dest='yaml',help='LinkML yaml file')
    parser.add_argument('--sha',dest='sha',help='SHA256 hash for file')
    args = parser.parse_args()
    

    fspath = args.figshare
    with open(args.yaml,'r') as file: data = yaml.safe_load(file)
    records = make_recordlist(data)
    distribution = make_filelist(args.figshare,args.sha)
    croissant = make_croissant(records,distribution)
    outfile = open('coderdata.json','w')
    json.dump(croissant, outfile, indent=4,ensure_ascii=True)
    outfile.close()

def make_recordlist(data):
    '''
    makes record list for croissant metadata object
    '''
    reclist=[] ##this is essentially the record set

    slot_desc=data['slots'] ##dictionary for slot info
    enums = data['enums']
    
    ##frst, every class should
    for tab in data['classes'].keys():
      #  if tab !='Sample':
      #      continue
        
        fields=[] ##every table definition gets its own fields
        keys = [] ##keys refer to fields

        if tab=='Gene': ###we need to treat this differently
            source='fileObject'
        else:
            source='fileSet'

        cr_id = tab
        if tab=='Copy Number':
            cr_id='CopyNumber'
        if tab=='Drug Descriptor':
            cr_id='DrugDescriptor'

        for slo in data['classes'][tab]['slots']:
            ##need to look up metadata for slot to format properly
            keys.append({'@id':tab+'/'+slo})
            desc=slot_desc[slo]['description']
            dtype='sc:Text'
            if 'range' in slot_desc[slo].keys():
                if slot_desc[slo]['range'] =='float':
                    dtype='sc:Float'
                elif slot_desc[slo]['range'] =='integer':
                    dtype='sc:Integer'
            nd = {
                '@type':'cr:Field',
                '@id': cr_id+'/'+slo,
                'description': desc,
                'dataType': dtype,
                'source': {source: {'@id': cr_id},
                           'extract': {'column': slo}
                           }
                }
            fields.append(nd)
        for attr in data['classes'][tab]['attributes'].keys():
            print(tab+'/'+attr)
            keys.append({'@id':tab+'/'+attr})
            ##attribute metdata should be inside
            aobj=data['classes'][tab]['attributes'][attr]
            desc=aobj['description']
            dtype='sc:Text'                     
            if 'range' in aobj.keys():
                if aobj['range'] =='float':
                    dtype='sc:Float'
                elif aobj['range'] =='integer':
                    dtype='sc:Integer'

            nd = {
                '@type':'cr:Field',
                '@id': cr_id+'/'+attr,
                'description': desc,
                'dataType': dtype,
                'source': {source: {'@id': cr_id},
                           'extract': {'column': attr}
                           }
                }
            fields.append(nd)
        reclist.append({
            '@type':'cr:RecordSet',
            '@id':tab,
            'key': keys,
            'field':fields,
            })
    return reclist



def make_filelist(figshare,sha):
    '''
    makes file list for distribution part of croissant metadata file
    '''
    dist = [
	    {
	        "@type": "cr:FileObject",
	        "@id": "coderdata-zip",
	        "contentUrl":figshare,
	        "encodingFormat":"application/zip",
	        "sha256": sha
	    },
	    {
	        "@type": "cr:FileObject",
	        "@id": "Gene",
	        "containedIn": {"@id": "coderdata-zip"},
	        "contentUrl" : "./genes.csv",
	        "description": "Gene identifier file",
	        "encodingFormat": "text/csv",
                "sha256": ""
	    },
	    {
	        "@type": "cr:FileSet",
	        "@id": "Sample",
	        "description": "Sample identifier files for each dataset",
	        "containedIn": {"@id": "coderdata-zip"},
	        "encodingFormat": "text/csv",
	        "includes": "*samples.csv"
	    },
	    {
	        "@type": "cr:FileSet",
	        "@id":"Drug",
	        "description":"Drug identifiers and SMILES string",
	        "containedIn": {"@id": "coderdata-zip"},
	        "encodingFormat": "text/tsv",
	        "includes": "*drugs.tsv"
	    },
	    {
	        "@type": "cr:FileSet",
	        "@id": "Transcriptomics",
	        "description": "Transcriptomics files",
	        "containedIn": {"@id": "coderdata-zip"},
	        "encodingFormat": "text/csv",
	        "includes": "*transcriptomics.csv.gz"
	    },
	    {
	        "@type": "cr:FileSet",
	        "@id": "Proteomics",
	        "description": "Proteomics files",
	        "containedIn": {"@id": "coderdata-zip"},
	        "encodingFormat": "text/csv",
	        "includes": "*proteomics.csv.gz"
	    },
	    {
	        "@type": "cr:FileSet",
	        "@id": "CopyNumber",
	        "description": "Copy number datafiles",
	        "containedIn": {"@id": "coderdata-zip"},
	        "encodingFormat": "text/csv",
	        "includes": "*copy_number.csv.gz"
	    },
	    {
	        "@type": "cr:FileSet",
	        "@id": "Mutations",
	        "description": "Mutation files",
	        "containedIn": {"@id": "coderdata-zip"},
	        "encodingFormat": "text/csv",
	        "includes": "*mutations.csv.gz"
	    }
        ]
    return dist

def make_croissant(recordlist,distribution):
    '''
    provides the context and mapping for the croissant header
    '''
    croissant={
        "@context": {
	    "@language": "en",
	    "@vocab": "https://schema.org/",
	    "citeAs": "cr:citeAs",
	    "column": "cr:column",
	    "conformsTo": "dct:conformsTo",
            "containedIn": "cr:containedIn",
            "description":"cr:description",
	    "cr": "http://mlcommons.org/croissant/",
	    "rai": "http://mlcommons.org/croissant/RAI/",
	    "data": {
	        "@id": "cr:data",
	        "@type": "@json"
	    },
	    "dataType": {
	        "@id": "cr:dataType",
	        "@type": "@vocab"
	    },
            "dct": "http://purl.org/dc/terms/",
	    "examples": {
	        "@id": "cr:examples",
	        "@type": "@json"
	    },
	    "extract": "cr:extract",
	    "field": "cr:field",
	    "fileProperty": "cr:fileProperty",
	    "fileObject": "cr:FileObject",
	    "fileSet": "cr:FileSet",
	    "includes": "cr:includes",
	    "isLiveDataset": "cr:isLiveDataset",
	    "jsonPath": "cr:jsonPath",
	    "key": "cr:key",
	    "md5": "cr:md5",
	    "format": "cr:format",
	    "ml": "http://mlcommons.org/schema/",
	    "parentField": "cr:parentField",
	    "path": "cr:path",
	    "recordSet": "cr:recordSet",
	    "references": "cr:references",
	    "regex": "cr:regex",
	    "repeated": "cr:repeated",
	    "replace": "cr:replace",
	    "sc": "https://schema.org/",
	    "separator": "cr:separator",
	    "source": "cr:source",
	    "subField": "cr:subField",
	    "transform": "cr:transform"
        },
        "@type": "sc:Dataset",
        "name": "coderdata",
        "conformsTo": "http://mlcommons.org/croissant/1.0",
        "description": "CoderData is a benchmark dataset designed for AI/ML modeling of drug sensitivity via omics measurements.",
        "citeAs": "@Article{jac24, author = \"Jeremy Jacobson, Chang In Moon, Sydney Schwartz, Belinda Garana, Ryan Weil, and Sara Gosline}\"",
        "license": "BSD-3",
        "url": "https://pnnl-compbio.github.io/coderdata/",
        "distribution": distribution,
        "recordSet": recordlist        
    }
    return croissant
    
    
if __name__=='__main__':
    main()
