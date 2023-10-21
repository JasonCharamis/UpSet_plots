import os
import re
import argparse
import pandas as pd
from upsetplot import UpSet
from matplotlib import pyplot
from itertools import combinations


def open_DE_file(filename):
    geneids = {}
    
    with open(filename, "r") as file:
        lines = file.readlines()

        for line in lines:
            if not re.search("Geneid|regulated", line):
                cols = line.split('\t')
                geneids[cols[1]] = line

    return geneids


def generate_lists_of_DE_genes(filenames):
    filename_regex = re.compile("../|.P1e-3_C2.DE.annotated.plus_orthology.sorted.tsv")

    subsets = []

    for filename in filenames:
        geneids = open_DE_file(filename)
        file_geneids_dict = {}

        file_geneids_dict[re.sub(filename_regex, "", filename) + "_upregulated"] = []
        file_geneids_dict[re.sub(filename_regex, "", filename) + "_downregulated"] = []

        for geneid in geneids.keys():
            if float(geneids[geneid].split('\t')[2]) > 0:
                file_geneids_dict[re.sub(filename_regex, "", filename) + "_upregulated"].append(geneid)
            else:
                file_geneids_dict[re.sub(filename_regex, "", filename) + "_downregulated"].append(geneid)

        for comparisons, genes in file_geneids_dict.items():
            if genes:
                subsets.append(file_geneids_dict)
            else:
                print(comparisons + " is empty!")

    return subsets
                
              
def upset_plots(input_files, outfile="UpSet_plot.png", DE=True):
    sets = []
    
    if DE:
        input_dict = {name: geneids for input_dict in generate_lists_of_DE_genes(input_files) for name, geneids in input_dict.items()}

        print(input_dict)

        if input_dict:
            names = [re.sub("_geneids", "", key) for key in input_dict.keys()]
            sets = list(input_dict.values())
        
        else:
            print("Problem with input dictionaries.")
            
    else:
        names = [re.sub(r"\.P1e-3_C2.DE.annotated.plus_orthology.sorted.*", "", os.path.splitext(os.path.basename(file_path))[0]) for file_path in input_files]
        sets = input_files
        
    all_elems = list(set().union(*sets))  # Unpack sets, find the unique elements and save them into a list
    df = pd.DataFrame([[e in st for st in sets] for e in all_elems], columns=names)  # Check if each of the unique elements is found in each subset
    df_counts = df.groupby(names).size()  # Group based on presence and compute size
    UpSet(df_counts, orientation='horizontal', subset_size='sum', show_counts=True).plot()
    pyplot.savefig(outfile)



def kargs():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Create Upset plots from a collection of files.')
    parser.add_argument('-DE, --DE_files', type=str, help='DE file to generate subsets of commonly up-regulated, commonly down-regulated and other interesting categories of genes.')
    parser.add_argument('-ls, --lists_of_strings', type=str, help='Provided lists of strings, in the case of not provided DE files.')
    parser.add_argument('-out, --outfile', type=str, help='Output png file to print UpSet plot.')
    
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()


def main():
    
    parser = argparse.ArgumentParser(description='Script to construct UpSet plots from DE files or lists of strings.')
    args = kargs()

    if not args.DE_files and not args.lists_of_strings:
        print ("Please provide a file as input.")
           
    elif args.DE_files and not args.lists_of_strings:
        if args.outfile:
          upset_plots ( args.DE_files, DE = True, outfile = args.outfile )
            
    elif not args.DE_files and args.lists_of_strings:
        if args.outfile:
          upset_plots ( args.lists_of_strings, DE = False, outfile = args.outfile )
        
            
if __name__ == "__main__":
    main()
