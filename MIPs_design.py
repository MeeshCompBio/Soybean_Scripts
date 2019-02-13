"""This sctip will take a gff file as input and parce out the longest
   transcripts  in gff format for formatting for exon mapping"""

import sys
import getopt
from collections import OrderedDict


def usage():
    print("""\n
        This is the usage function:
        python3 Parse_GFF_LongestTranscript.py -g <gff_file> -o <output_file>
            -g or gff_file    :The name of the gff file you want to use
            -o or output_file :The name of your GO network object
            -t or --type      :Scan the gff for low or high priority (default)
                               genes. If this flag is used then it will search
                               for low priority genes using raw gff file. High
                               priority search requires a parsed gff of only
                               gene and exon fields for whatever transcript
                               representation you choose. Check out my github
                               repo https://github.com/MeeshCompBio/Basic_MSI_
                               Utilities for a script to pull out the longest
                               representative transcript that can be used on
                               a raw gff.
        \n""")


try:
    opts, args = getopt.getopt(sys.argv[1:], "g:o:th", ["gff_file=",
                                                        "output_file=",
                                                        "type="
                                                        "help"])
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)
Type = "HP"
for opt, arg in opts:
    if opt in ("-g", "--gff_file"):
        gff_file = arg
    elif opt in ("-o", "--output_file"):
        output_file = arg
    elif opt in ("-t", "--type"):
        Type = "LP"
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"

# open file for reading, writing and ID len
in_file = open(gff_file, "r")
out_file = open(output_file, "w")


# maybe include type for high profile or not
# print the header
print("MIP_id",
      "MIP_Chr",
      "MIP_start",
      "MIP_stop",
      "Glyma",
      "Glyma_start",
      "Glyma_stop",
      "Glyma_size",
      "Exon_name",
      "Exon_start",
      "Exon_stop",
      "exon_length",
      sep="\t")


# What information do I want to store
# MIPs_ID, Gene, Chrom, start, stop, type
mip_id = 0
if Type == "HP":
    # initialize the first round  of dicts even though it will get looped over
    for row in in_file:
        line = row.rstrip()
        # Keep all of the header files
        if line.startswith("#"):
            continue
        else:
            # We are just looking for information then printing the entire row
            # so we are not going to strip newline
            fields = line.split("\t")
            # If we find a gene, this signifies the step to find the longest
            # transcript for the prevoius gene
            if fields[2] == "gene":
                # This is to prevent errors on first loop with
                gene = fields[8].split(";")[0].split("=")[1]
                gene_chr = fields[0]
                gene_start = fields[3]
                gene_stop = fields[4]
                gene_length = int(gene_stop)-int(gene_start)
                exon_num = 0
            elif fields[2] == "exon":
                exon_num += 1
                exon_name = fields[8].split(";")[0].split("=")[1]
                mip_id += 1
                exon_chr = fields[0]
                exon_start = fields[3]
                exon_stop = fields[4]
                exon_length = int(exon_stop) - int(exon_start)
                half_dist = int(exon_length/2)

                # not the best way to do this I should refactor later
                if exon_length <= 50:
                    # Just take the midpoint
                    print(mip_id,
                          gene_chr,
                          ((int(exon_start)+half_dist)-100),
                          ((int(exon_start) + half_dist) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")

                elif exon_length <= 150:
                    # take the intron and exon boundaries
                    print(mip_id,
                          gene_chr,
                          (int(exon_start)-100),
                          (int(exon_start) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")
                    mip_id += 1
                    print(mip_id,
                          gene_chr,
                          (int(exon_stop)-100),
                          (int(exon_stop) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")

                elif exon_length <= 300:
                    # Take intron/exon bondaries, plus midpoint
                    print(mip_id,
                          gene_chr,
                          (int(exon_start)-100),
                          (int(exon_start) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")
                    mip_id += 1
                    print(mip_id,
                          gene_chr,
                          ((int(exon_start)+half_dist)-100),
                          ((int(exon_start) + half_dist) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")
                    mip_id += 1
                    print(mip_id,
                          gene_chr,
                          (int(exon_stop)-100),
                          (int(exon_stop) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")

                else:
                    # Take intron/exon boundaries then tiled every 150 bp
                    mip_exon_idx = 0
                    print(mip_id,
                          gene_chr,
                          (int(exon_start)-100),
                          (int(exon_start) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")
                    mip_id += 1

                    while (mip_exon_idx+150) <= exon_length:
                        mip_exon_idx += 150
                        print(mip_id,
                              gene_chr,
                              ((int(exon_start)+mip_exon_idx)-100),
                              ((int(exon_start) + mip_exon_idx) + 100),
                              gene,
                              gene_start,
                              gene_stop,
                              gene_length,
                              exon_name,
                              exon_start,
                              exon_stop,
                              exon_length,
                              sep="\t")
                        mip_id += 1
                    print(mip_id,
                          gene_chr,
                          (int(exon_stop)-100),
                          (int(exon_stop) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")

elif Type == "LP":
    def output_exon(list_of_lists):
        """
            Take a list of lists of exon candidates and see which one
            is larger. We only compare the first two in the list of lists so
            it assumes that it is sorted beforehand
        """
        if len(list_of_lists) == 0:
            return("Exons don't overlap")
        elif len(list_of_lists) == 1:
            exon_1 = list_of_lists[0]
            exon_name = exon_1[0]
            exon_start = exon_1[1]
            exon_stop = exon_1[2]
            exon_length = exon_1[3]
            half_dist = int(exon_length/2)
        else:
            exon_1 = list_of_lists[0]
            exon_2 = list_of_lists[1]

            if int(exon_1[3]) >= int(exon_2[3]):
                exon_name = exon_1[0]
                exon_start = exon_1[1]
                exon_stop = exon_1[2]
                exon_length = exon_1[3]
                half_dist = int(exon_length/2)
            else:
                exon_name = exon_2[0]
                exon_start = exon_2[1]
                exon_stop = exon_2[2]
                exon_length = exon_2[3]
                half_dist = int(exon_length/2)
        return(exon_name,
               exon_start,
               exon_stop,
               exon_length,
               half_dist,
               )

    # make functions unique to this path
    def multi_match_exons(ordered_dict):
        """take in an ordered dictionary and make sure the first two exons
           are represented in at least 50% of transcripts, if not add the next
           exon in line (i.e. 3rd) until that result is satisfied. Determine
           which of the two exon is larger and return the metadata for that
           exon
        """
        size = len(ordered_dict)
        transcript_1 = ordered_dict.popitem(last=False)[1]
        transcript_2 = ordered_dict.popitem(last=False)[1]
        multi_exon = []
        new_size = size-2
        # check if there is only two exons in TS
        if new_size == 0:
            # check to see which exons are in both TS
            for item in range(len(transcript_1)):
                # see if exon size matches across transcripts
                if transcript_1[item][1:] in [x[1:] for x in transcript_2]:
                    multi_exon.append(transcript_1[item])

            return(output_exon(multi_exon))

        else:
            multi_transcript = [transcript_1, transcript_2]
            while len(ordered_dict) != 0:
                multi_transcript.append(ordered_dict.popitem(last=False)[1])

            exon_sizes = []
            for ts in multi_transcript:
                for exon in ts:
                    if exon[1:] not in exon_sizes:
                        exon_sizes.append(exon[1:])

            exon_in_each_ts = []
            for ex in exon_sizes:
                counter = 0
                for ts in multi_transcript:
                    if ex in [x[1:] for x in ts]:
                        counter += 1
                        candidate_exon = ts[[x[1:] for x in ts].index(ex)]
                if counter == size:
                    exon_in_each_ts.append(candidate_exon)
            # Sort them by exon numer at the end of string
            exon_in_each_ts = sorted(exon_in_each_ts,
                                     key=lambda x: int(x[0].split(".")[-1]))

            return(output_exon(exon_in_each_ts))

    mip_id += 1
    Transcript_Compare = OrderedDict()
    Gff_line_exon = dict()
    for row in in_file:
        line = row.rstrip()
        # Keep all of the header files
        if line.startswith("#"):
            continue
        else:
            # We are just looking for information then printing the entire row
            # so we are not going to strip newline
            fields = line.split("\t")
            # If we find a gene, this signifies the step to find the longest
            # transcript for the prevoius gene
            if fields[2] == "gene":
                # This is to prevent errors on first loop with
                if bool(Transcript_Compare):
                    # Return key with largest value
                    if len(Transcript_Compare) == 1:
                        transcript = Transcript_Compare.popitem(last=False)[1]
                        if len(transcript) == 1:
                            exon_name = transcript[0][0]
                            exon_start = transcript[0][1]
                            exon_stop = transcript[0][2]
                            exon_length = transcript[0][3]
                            half_dist = int(exon_length/2)
                        else:
                            exon_1 = transcript[0]
                            exon_2 = transcript[1]
                            if int(exon_1[3]) >= int(exon_2[3]):
                                exon_name = exon_1[0]
                                exon_start = exon_1[1]
                                exon_stop = exon_1[2]
                                exon_length = exon_1[3]
                                half_dist = int(exon_length/2)
                            else:
                                exon_name = exon_2[0]
                                exon_start = exon_2[1]
                                exon_stop = exon_2[2]
                                exon_length = exon_2[3]
                                half_dist = int(exon_length/2)
                    else:
                        multi_exon = multi_match_exons(Transcript_Compare)
                        if "Exons don't overlap" in multi_exon:
                            print("mip_id",
                                  "Exons don't overlap or " +
                                  "exons not present in each transcripts",
                                  "XXXXX",
                                  "XXXXX",
                                  gene,
                                  gene_start,
                                  gene_stop,
                                  gene_length,
                                  "XXXXX",
                                  "XXXXX",
                                  "XXXXX",
                                  "XXXXX",
                                  sep="\t")
                        else:
                            exon_name = multi_exon[0]
                            exon_start = multi_exon[1]
                            exon_stop = multi_exon[2]
                            exon_length = multi_exon[3]
                            half_dist = multi_exon[4]

                    # print to the output depending on the options
                    print(mip_id,
                          gene_chr,
                          ((int(exon_start)+half_dist)-100),
                          ((int(exon_start) + half_dist) + 100),
                          gene,
                          gene_start,
                          gene_stop,
                          gene_length,
                          exon_name,
                          exon_start,
                          exon_stop,
                          exon_length,
                          sep="\t")
                    mip_id += 1

                else:
                    gene = fields[8].split(";")[0].split("=")[1]
                    gene_chr = fields[0]
                    gene_start = fields[3]
                    gene_stop = fields[4]
                    gene_length = int(gene_stop)-int(gene_start)
    ##########################################################################
                    # # This is just for testing, uncomment if you want to see
                    # # the length outputs for each transcript and the longest
                    # for k, v in Transcript_Compare.items():
                    #     print(k, v)
                    # print("Longest=", Longest_transcript)
    ##########################################################################

                # reinitialize/empty dicts for next gene
                Transcript_Compare = OrderedDict()
                Gff_line_exon = dict()
                # grab information for the next set of genes
                gene = fields[8].split(";")[0].split("=")[1]
                gene_chr = fields[0]
                gene_start = fields[3]
                gene_stop = fields[4]
                gene_length = int(gene_stop)-int(gene_start)

            elif fields[2] == "exon":
                # Get the size of the exon
                exon_start = int(fields[3])
                exon_stop = int(fields[4])
                exon_size = exon_stop - exon_start
                # Pull out and split information from ID category
                exon_ID = fields[8].split(";")[0].split("=")[1]
                # Pull out the transcript ID
                Transcript = exon_ID.split(".exon")[0]
                Parent = fields[8].split(";")[1].split("=")[1]

                # Look to see if that transcript is already in the dict
                if Transcript in Transcript_Compare.keys():
                    Transcript_Compare[Transcript].append([exon_ID,
                                                           exon_start,
                                                           exon_stop,
                                                           exon_size
                                                           ]
                                                          )
                    # Store the raw row information using the same key
                    Gff_line_exon[Transcript].append(row)
                else:
                    # if not make a new one
                    Transcript_Compare[Transcript] = [[exon_ID,
                                                       exon_start,
                                                       exon_stop,
                                                       exon_size
                                                       ]
                                                      ]
                    Gff_line_exon[Transcript] = [row]
    # This is for the last gene/exon information that is stored
    if bool(Transcript_Compare):
        # Return key with largest value
        if len(Transcript_Compare) == 1:
            transcript = Transcript_Compare.popitem(last=False)[1]
            # print(len(transcript))
            if len(transcript) == 1:
                exon_name = transcript[0][0]
                exon_start = transcript[0][1]
                exon_stop = transcript[0][2]
                exon_length = transcript[0][3]
                half_dist = int(exon_length/2)
            else:
                exon_1 = transcript[0]
                exon_2 = transcript[1]
                if int(exon_1[3]) >= int(exon_2[3]):
                    exon_name = exon_1[0]
                    exon_start = exon_1[1]
                    exon_stop = exon_1[2]
                    exon_length = exon_1[3]
                    half_dist = int(exon_length/2)
                else:
                    exon_name = exon_2[0]
                    exon_start = exon_2[1]
                    exon_stop = exon_2[2]
                    exon_length = exon_2[3]
                    half_dist = int(exon_length/2)
        else:
            multi_exon = multi_match_exons(Transcript_Compare)
            exon_name = multi_exon[0]
            exon_start = multi_exon[1]
            exon_stop = multi_exon[2]
            exon_length = multi_exon[3]
            half_dist = multi_exon[4]
        # print to the output depending on the options
        print(mip_id,
              gene_chr,
              ((int(exon_start)+half_dist)-100),
              ((int(exon_start) + half_dist) + 100),
              gene,
              gene_start,
              gene_stop,
              gene_length,
              exon_name,
              exon_start,
              exon_stop,
              exon_length,
              sep="\t")
        mip_id += 1

# close files for safety
in_file.close()
out_file.close()
