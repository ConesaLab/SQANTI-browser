#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Alerjandro Paniagua
# Script to add the SQANTI3 cantegories color, the structural category and the
# associated transcript

import os
import csv


def get_RGB(structural_category: str)-> str:
    '''
    Function to get the RGB color corresponding to a structutal category
    '''
    # Structural category corrsponding colors
    cat_palette = {"full-splice_match":"6BAED6", 
                   "incomplete-splice_match":"FC8D59", 
                   "novel_in_catalog":"78C679", 
                   "novel_not_in_catalog":"EE6A50", 
                   "genic":"969696",
                   "antisense":"66C2A4",
                   "fusion":"DAA520",
                   "intergenic" : "E9967A", 
                   "genic_intron":"41B6C4"}
    
    # Get HEX code
    hex = cat_palette[structural_category]
    
    # Covert HEX to RGB (it's what Genome Browser uses)
    rgb = []
    for i in (0, 2, 4):
        decimal = int(hex[i:i+2], 16)
        rgb.append(str(decimal))

    return ",".join(rgb)


def parse_class(class_file: str)-> dict:
    '''
    Function to parse the stucural category and associated transcript into a
    dictionary

    Args:
        class_file (str): path to the classification file

    Output:
        class_dict (dict): dictionary with the relevant information for each
            transcript
    '''
    class_dict = dict()
    with open(class_file, "r") as f_in:
        csv_reader = csv.DictReader(f_in, delimiter="\t")
        class_dict = dict()
        # For each transcript save the SC, associated gene and get RGB color
        for line in csv_reader:
            class_dict[line["isoform"]] = (line["structural_category"], 
                                           line["associated_transcript"], 
                                           get_RGB(line["structural_category"]))

    return class_dict
            

def main():
    transcript_id_bed = 3
    rgb_field = 8

    bed_file, class_file = snakemake.input
    bed_out = snakemake.output[0]
    # Get the transcritpts info
    class_dict = parse_class(class_file)

    # Get the bed file name wo extension
    bed_name, extension = os.path.splitext(bed_file)


    with open(bed_file, "r") as f_in, open(bed_out, "w") as f_out:
        csv_reader = csv.reader(f_in, delimiter="\t")
        csv_writer = csv.writer(f_out, delimiter="\t")

        for line in csv_reader:
            transcript_id = line[transcript_id_bed]
            if transcript_id in class_dict.keys():
                structural_category, associated_transcript, RGB = class_dict[transcript_id]
                line[rgb_field] = RGB
                line.append(structural_category)
                line.append(associated_transcript)
            else:
                line[rgb_field] = "255,255,255"
                line.append("rescued")
                line.append(transcript_id)              
            csv_writer.writerow(line)

if __name__ == "__main__":
    main()