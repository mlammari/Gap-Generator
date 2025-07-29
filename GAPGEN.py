# GAPGEN - edits 

# Sample Command: python Gap_Generator.py -in /Users/meenaammari/Documents/Paleogenomics_Lab/seq/ls_orchid.fasta -ng 3 -out new_orchid.fasta -s orchid

# This program develops a script that introduces gaps to a reference genome:
# - Introduce parameters
# - Add number of gaps and size of gaps
# - No gaps at the beginning and end of a sequence 
# - Keep track of the sequences I am substiuting for gaps 
# - Note: gaps = sequences of N's that have a specific size

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import chain
from datetime import date
import random, argparse
import numpy as np 
import matplotlib.pyplot as plt
import sys, os
import psutil

#def parse_data(file):
#    """ This function takes the input fasta file and condenses it down to the entire genome by removing the headers and creating a list with one element. That 
#        element being the entire genome. 

 #   Args:
  #      file (_file_): The input fasta file. 

  #  Returns:
  #      _type_: _description_
    #parsed_data = SeqIO.parse(file, "fasta")
#    gen_list = []
#    new_gen_list = []
#    counter = 0
#    with open(file, "r") as gen:
#        for line in gen:
#            if line[0] == ">":
#                continue
#            else:
#                gen_list.append(line)
##        for s in gen_list:
 #           if s == "\n":
    #            gen_list.remove(s)
    #            counter += 1
    #        else:
    #            new_gen_list.append(gen_list[counter].strip())
            #    #for letter in string:
            #    #    if letter == "\n"
            #    #        letter
    #        #counter += 1
    #dna = "".join(new_gen_list)
    
    #return dna 
    
    
    #comp_list = "".join(gen_list)
    #print(comp_list)

def computeLength(file_input): # will also parse through data for testing purposed
    data = SeqIO.parse(file_input, "fasta")
    rec_count = 0
    rec_list = []
    for rec in data:
        rec_count += 1
        rec_list.append(rec)

    #for record in rec_list:
    #print(rec_list[0].seq)
        
    
    return rec_count, rec_list

                
def percentiile(file_input, threshold):
    """ This function does two main tasks. the first task is to compute the length at the user-defined percentile of all the data. This length will be 
        returned to main to determine whether a gap should be inserted in each sequence. The second task is to create a histogram showing the read distribution
        of all the reads in the genome; x-axis = length of each sequence, y-axis = number of reads.

    Args:
        file_input (str): The fasta file that the user input. 
        threshold (int): The percentile threshold the user wants to set.

    Returns:
        int : This function returns the number at the percentile that user defined. 
    """
    information = SeqIO.parse(file_input, "fasta")
    lengths = []

    [lengths.append(len(record.seq)) for record in information]

    lens_array = np.array(lengths)

    num_percentile = np.percentile(lens_array, threshold)

    plt.hist(lens_array, bins = 10, edgecolor = "black")
    plt.title("Sequence Distribution")
    plt.xlabel("Sequence Length")
    plt.ylabel("Number of Sequences")
    plt.savefig("percentile.jpeg")

    return num_percentile
    

def sequenceLength(sequence):
    """ This is a parameter that ensures that the input sequence is less that 10000000  bases long  
    Args:
        sequence (_Bio.Seq.Seq_): This is the input sequence that we will conduct our check on 

    Raises:
        ValueError: If the sequence is greater than 10000000 bases long, the funtion raises a Value Error 

    Returns
        _boolean_: If the function returns true, then the input sequence is less than 100 bases long 
    """
    if len(sequence) > 100000000:
        return False
        #raise ValueError\
        #("The length of your sequence is too long")
    else:
        return True 
def GC_Content(sequence):  
    """ This is a parameter that checks to see if the GC fraction of the sequence is less than 55%
    Args:
        sequence (_Bio.Seq.Seq_): This is the input sequence that we will conduct our check on

    Raises:
        ValueError: If the sequence has a GC fraction of over 55%, the function raises a Value Error

    Returns:
        _boolean_: If the function returns true, then the input sequence has a GC fraction that is less than 55 %
    """
    my_seq = Seq(sequence)
    if gc_fraction(my_seq) > 0.55:
        return False
    else: 
        return True
    
def assignGaps(file_input, num_gaps):
    parse_data = SeqIO.parse(file_input, "fasta")
    records = []
    gap_rec = 0

    for rec in parse_data:
        records.append(rec.id)
    
    if num_gaps >= len(records):
        gap_rec = num_gaps // len(records)
                    #100        #for orchid == 94
    else:
        gap_rec = 1
    
    remainder = num_gaps % len(records)

   
    return gap_rec, remainder

def snip_seq(sequence, gap_size): 
    """ This function takes in the sequence and clips off the begininng/end of those sequences while also saving those snipped off ends to be 
        concatenated later. The length of the snipped ends will be equal to the gap size. 

    Args:
        sequence (_Bio.Seq.Seq_): This is the input sequence that we will clip the ends off of
        gap_size (_int_): The size of the gap, also the length of the sequence that will be snipped off

    Returns:
        beg_seq, end_seq, updated_seq : This function returns three variables:
        
        - beg_seq: the beginning sequence that we snipped off
        - end_seq: the end sequence that we snipped off
        - updated_seq: the sequence we will actually interact with. In other words, the sequence will implement gaps in

    """
    seq_list = sequence.split()

    str_seq = seq_list[0]
    new_seq = str_seq.split()
    
    updated_seq = seq_list[0][gap_size:(len(sequence)-gap_size)] 
    beg_seq = []            
    end_seq = []

    for base in range(gap_size):            
        beg_seq.append(seq_list[0][base])   
    
    for base in range(gap_size):    
        end_seq.append(seq_list[0][-base])  
    
    return beg_seq, end_seq, updated_seq    

def getNumGaps(num_gaps, ng_fixed):
    """ This function produces the number of gaps that will be produced in the sequence. Whether or not the number of gaps will be fixed or not will depend on the user's
        specifications. 

    Args:
        num_gaps (_int_): The highest number of gaps that can be produced in the sequence. 
        ng_fixed (_bool_): Determines whether or not the number of gaps will be fixed or vary-able across the sequences. 

    Returns:
        _int_: The number of gaps that will be placed into the sequence. 
    """
    if ng_fixed:
        poss_gaps = random.randint(1,num_gaps)
    else:
        poss_gaps = num_gaps
    
    return poss_gaps

def getSequenceSize(seq_size, seq_fixed):
    """ This function produces the size of the sequence that will be removed from the sequence. Whether or not the number of gaps will be fixed or not will depend on 
    the user's specifications. 

    Args:
        seq_size (_int_): The maximum size of the sequence that will be removed from the orginal contig/scaffold. 
        seq_fixed (_bool_): Determines whether or not the size of the sequence will be fixed or vary-able acrros the sequence. 

    Returns:
        _int_: The size of the sequence that will be replaced with a gap. 
    """

    if seq_fixed: 
        start_rand = random.randint(0,seq_size)
    else:
        start_rand = seq_size

    return start_rand

def validCoordinate(updated_seq, poss_gaps, start_rand): 
    """ This function generates a list of valid coordinates of the sequences that will be removed from the orignal sequence and replaced with gaps.
        The parameters for this function are the updated sequence, the highest number of gaps that could or would be placed (depending on user specifications)
        and the actual size of the sequence we will replace (again based on user specifications)

    Args:
        updated_seq (_str_): The updated sequence after we have snipped off the ends of the sequence
        poss_gaps (_int_): The highest number of gaps that could or would be placed (depending on user specifications)
        start_rand (_int_): The actual size of the sequence we will replace (again based on user specifications)

    Returns:
        _list_: gap_coordinates - a list of of the coordinates that will be replaced with a gap
    """
    gap_coordinates = []

    start_coordinate = random.randint(0,len(updated_seq)-start_rand)
    gap_coordinates.append([start_coordinate, start_coordinate + start_rand])

    for x in range(1, poss_gaps): # Thus far I have kept the range for this for loop to go from 1 to poss_gaps. This may change as the code develops. 
        isValid = False
        while isValid == False:
            coordinate = random.randint(0,len(updated_seq) - start_rand) #rand_size
            count = 0
            for tup in gap_coordinates:
                lower_bound = tup[0]
                upper_bound = tup[1]                                                                #rand_size
                if (coordinate not in range(lower_bound -1, upper_bound + 1)) and ((coordinate + start_rand) not in range(lower_bound -1, upper_bound+1)):
                    count += 1
                    if lower_bound in range(coordinate, coordinate + start_rand): #rand_size
                        count -= 1 
                        break
                if count == len(gap_coordinates):  
                    gap_coordinates.append([coordinate, coordinate + start_rand]) #rand_size
                    isValid = True
               
    return gap_coordinates
    

def create_gaps(sequence, gap_coordinates, gap_size): 
    """ This function takes in the updated sequence and actually generates the gaps at the random coordinates we generated in the 
        validCoordinate function. 

    Args:
        sequence (_Bio.Seq.Seq_): The updated sequence (the ends snipped off) we will implement gaps in 
        coordinates (_list_): A list containing tuples of the random coordinates at which the gaps will be placed
        gap_size (_int_): The size of the gap. The number of N's that will be placed in the new DNA sequenced. 
    """  
    gap_coordinates.sort()
    

    seq_track = []
    new_seq = list(sequence)

    for i in range(len(gap_coordinates)): 
        first_gap = sequence[gap_coordinates[i][0]:gap_coordinates[i][1]] 
        seq_track.append(first_gap)
        
    for i in range(len(gap_coordinates)): 
        start = gap_coordinates[i][0]
        end = gap_coordinates[i][1]
        gap = gap_size * "N"
        new_seq[start:end] = gap

    #print(seq_track)

    return new_seq, seq_track

def create_agp(file_input, agp_file, species):
    """ This function takes in the input file and the name of the agp file to generate the agp file that we will append the genome's information to

    Args:
        file_input (_str_): The path to the input FASTA file
        agp_file (_str_): The user-defined name of the agp file that will be created

    Returns:
        _str_: Returns the agp file
    """
    day = date.today()
    
    today = str(day)

    header = ["# agp_2.1\n",
              "#", file_input,"\n", 
              "# Assembly Name: ", species,"\n"
              "# Generated with Gap_Generator on ", today, "\n"]
               #maybe include command line prompt here

    with open(agp_file, "w") as agp:
        agp.write("".join(header))

    return agp_file

def append_agp(coor_dict, agp_file, seq_lens, gap_size, species):
    """ This function takes in the variables below and appends information to our agp file. The function iterates through each value in our 
    coor_dict dictionary, and appends either a sequence line or a gap line.  


    Args:
        coor_dict (_dict_): This is a dictionary that has the sequence id as keys, and its respecitve sequence coordinates as its values
        agp_file (_str_): This is the file we creates in create_agp(). In this function we take in this file and append our information to it
        seq_lens (_list_): contains a list of all the sequence lengths in our FASTA file
        gap_size (_int_): The gap size (sequence of N's) that will be used when writing the "New-Sequence" part of the agp file
        species (_str_): A user argument specifying the species name
    """
    counts = 1
    seq_counter = 1

    space = "  "

    with open(agp_file, "a") as agp:
        for key,seq_length in zip(coor_dict,seq_lens):
            if coor_dict[key] == "None":
                information = [species+"_"+str(counts)+"\t",
                               "1"+space+str(seq_length)+"\t",
                               str(seq_counter)+"\t",
                               "W\t",
                               key+"\t",
                               "1"+space+str(seq_length)+"\t",
                               "+",
                               "\n"]
                agp.write("".join(information))
                seq_counter = 1
                counts += 1
            else:
                new_list  = list(chain(*coor_dict[key]))  # I got this method from online: https://datagy.io/python-flatten-list-of-lists/
                #print(new_list)
                new_list.append(1)
                new_list.append(seq_length)

                new_list.sort()

                seq_coor = new_list

                beg = seq_coor[0]
                end = seq_coor[1]

                var_1 = end - 1
                
                initiator = [species+"_"+str(counts)+"\t",
                               str(beg)+space+str(var_1)+"\t",
                               str(seq_counter)+"\t",
                               "W\t",
                               key+"\t",
                               str(beg)+space+str(end-1)+"\t",
                               "+",
                               "\n"
                                ]
                
                agp.write("".join(initiator))

                seq_counter += 1

                for index in range(2, len(seq_coor), 2):
                    
                    pos_1 = seq_coor[index]
                    pos_2 = seq_coor[index+1]

                    new_pos = pos_2 - pos_1

                    if pos_1 > 1:
                        information = [species+"_"+str(counts)+"\t",
                               str(var_1 + 1)+space+str(var_1+1+gap_size)+"\t",
                               str(seq_counter)+"\t",
                               "N\t",
                               str(gap_size)+"  scaffold"+space+"no\t na",
                               "\n"]
                        var_1 += gap_size + 1
                        agp.write("".join(information))
                        seq_counter += 1
                    
                    information = [species+"_"+str(counts)+"\t",
                               str(var_1+1)+space+str(var_1+1+new_pos)+"\t",
                               str(seq_counter)+"\t",
                               "W\t",
                               key+"\t",
                               str(pos_1)+space+str(pos_2)+"\t",
                               "+",
                               "\n"
                                ]
                    seq_counter += 1
                    var_1 += 1 + new_pos
                    
                    
                    agp.write("".join(information))

                counts += 1

            seq_counter = 1
    return new_list

def create_bed(bed_file, coor_dict, identifier, flank):
    """ This file creates a BED file detailing all the locations of the gaps we implemented into our genome. It takes in three parameters described below. 

    Args:
        bed_file (_str_): The name of the BED file that we will create.
        coor_dict (_dict_): This is a dictionary that has the sequence id as keys, and its respecitve sequence coordinates as its values
        identifier (_str_): The identifier that will be used for each line in the BED file (will appear on the far right of each line in the file). 
    """
    iterate = 1
    with open(bed_file, "w") as bed:
        for record in coor_dict:
            if coor_dict[record] == "None":
                entry = [record+"\t"+"NA"+"\n"]
            else:
                #print(coor_dict[record])
                for gap in coor_dict[record]:
                    left = gap[0]-flank
                    right = gap[1]+flank
                    entry = [record+"\t"+str(left)+"\t"+str(right)+"\t"+identifier+str(iterate)+"\n"]
                    bed.write("".join(entry))
                    iterate += 1

def create_gap_fasta(gap_fasta, coor_dict, all_gaps, identifier):
    """ This function creates a FASTA file with all the sequences we replaced for gaps and its respective seq record. It takes in three parameters as explained below

    Args:
        gap_fasta (_str_): A file that the user names. This is the FASTA file that the function will write to 
        coor_dict (_dict_): This is a dictionary that has the sequence id as keys, and its respecitve sequence coordinates as its values
        all_gaps (_list_): This is a nested list containing all the sequences we replaced for gaps. The sequences in this list match the coordinates in coor_dict
        identifier (_str_): The identifier that will be used for each line in the BED file (will appear on the far right of each line in the file). 
    """
    item = 0
    counter = 1
    with open(gap_fasta, "w") as gapped:
        for record in coor_dict:
            while item < len(all_gaps):
                #print(all_gaps[item])
                for sequence in all_gaps[item]:
                    info = SeqRecord(seq=sequence, id=identifier+str(counter))
                    write = SeqIO.write([info], gapped, "fasta")
                    counter += 1
                item += 1


def flank_fasta(flank, flanks_fasta, sequence, id, gap_coordinates, mode, fcount):
    """ This function creates a fasta file (through each iteration as the main for loop proceeds through each record) containg both the sequence that was substituted
        for a gap, as well as the flanking sequences.

    Args:
        flank (_int_): The size of the flank. In this case the size of the flank will be less than or equal to the gap size. 
        flanks_fasta (_str_): The name of the fasta file that contains the flanked sequebces. 
        seq_input (_Bio.Seq.Seq_): The complete contig sequence.
        id (_strD_): The identifier for each sequence
        gap_coordinates (_list_): In this function, gap_coordinates is a list of the locations of all the gaps OF THE COMPLETE, UNEDITIED SEQUENCE. 
        mode (_str_): Determines which mode to write the fasta file to. 
    """ 
    
    with open(flanks_fasta, mode) as fd:
        for pair in gap_coordinates:
            seq = sequence[pair[0]-flank:pair[1]+flank]
            info = SeqRecord(seq=seq, id="flanks+"+id+str(fcount))
            writer = SeqIO.write([info], fd, "fasta")   
            fcount += 1
    return fcount

#def compile(file_input, threshold): 

def sampleNumberGaps(tota_gaps, num_sequences):
    """
    This function generates a dictionary with the number of gaps that will be generated per sequence.

    Args:
        total_gaps(int): total number of gaps that will be generated 
        sequences(list): list of sequences where gaps are going to be added
    Output:
        dictionary: indices are sequencesID, values are numner of gaps
    """
    #poss_gaps, remainder = assignGaps(file_input, num_gaps)
    gaps_list=[poss_gaps for recID in range(0,num_sequences)]
    for rem in range(0, remainder):
        gaps_list[rem]+=1
    return gaps_list

if __name__ == "__main__": 
    
    """
    This portion of the code performs all the necessary functions in order to do the following:
        - presents the arguments available to the user using argparse. These arguments include:
            - the path to the input file, the possible number of gaps the user wants to implement into each sequence (required)
            - the gap size (the number of N's the user will want to see if the final fasta file (required)
            - the path and the name of the new fasta file that will be generated (optional)
            - the probability that the algorithm will place a gap into the input sequence (optional)
            - the name of the agp file that will automatically be generated when the program is run
        - writes and appends the new, gapped (or not gapped) sequences to the new fasta file the user defined. 
        - calls the agp file creation function in order to generate the agp file

    """

    parser = argparse.ArgumentParser(description= "Generates Random Gaps")

    parser.add_argument("--input_file", "-in", type=str, help= "Path to file", required=True)

    parser.add_argument("--percentile", "-pr", metavar= "PERCENTILE", type=int, help= "The percentile this program will act on or above", default=0)

    parser.add_argument("--num_gaps", "-ng", metavar= "NUM_GAPS", type=int, help = "Highest possible number of gaps", required=True)

    parser.add_argument("--ng_fixed", "-fng", action="store_true", help="The number of gaps is set to a fixed value, provide the argument to allow the number of gaps to vary")
    
    parser.add_argument("--seq_size", "-sq", metavar= "ACTUAL_SIZE", type=int, help = "Actual size of the sequence that will be removed from the original reference", required=True)

    parser.add_argument("--seq_fixed", "-fsq", action="store_true", help ="The sequence size that will be removed is a fixed value, provide the argument to allow the sequence size to vary")  

    parser.add_argument("--gap_size", "-gs", metavar="GAP_SIZE", type = int, help="Gap size (number of N's) ", default=20)

    parser.add_argument("--flank_size", "-fs", metavar="FLANK_SIZE", type=int, help="The size of the flank; the sequence that preceeds and follows each gap)", default=20)

    parser.add_argument("--output_file", "-out", type= str, help= "The new fasta file you would like to generate ", required=True) # Might add default file name

    #parser.add_argument("--probability", "-p", metavar="PROBABILITY", type=int, help = "The percent chance the Gap Generator will generate a gap", default=100)

    parser.add_argument("--agp_file", "-a", type=str, help= "The name you would like for the agp file that will be generated", default="agp_file.agp")

    parser.add_argument("--prefix", "-s", metavar= "PREFIX", type=str, help="The prefix that will appear on the AGP file you create", required=True)

    parser.add_argument("--bed_file", "-b", type=str, help="The name of the bed file you would like to generate", default="bed_file.bed")

    parser.add_argument("--gap_fasta", "-gp", type=str, help="The name of the gap-fasta file you would like to generate", default="gap_fasta.fasta")

    parser.add_argument("--flanks_fasta", "-fsf", type=str, help="The name of the flanks fasta file that the program will generate", default="flanks_fasta.fasta")

    parser.add_argument("--gap_identifier", "-i", type=str, help="The name of the identifier that you would like to use to identfy the sequence in both the gap fasta file and the bed file", default="gap")

    args = parser.parse_args()
    
    file_input = args.input_file
    num_gaps = args.num_gaps
    ng_fixed = args.ng_fixed
    seq_size = args.seq_size
    seq_fixed = args.seq_fixed
    gap_size = args.gap_size
    #probability = args.probability
    agp_file = args.agp_file
    species = args.prefix
    bed_file = args.bed_file
    gap_fasta = args.gap_fasta
    identifier = args.gap_identifier
    flank = args.flank_size
    flanks_fasta = args.flanks_fasta
    threshold = args.percentile

    #parsed_data = SeqIO.parse(file_input, "fasta")
    num_percentile = percentiile(file_input, threshold)

    seq_lens = []
    start_ends = []
    all_gaps = []
    coor_dict = {}

    iterate = 1
    iter = 0
    counts = 0
    fcount = 1

    poss_gaps, remainder = assignGaps(file_input, num_gaps)
    total_records, record_list = computeLength(file_input)
    gaps_list = sampleNumberGaps(num_gaps, total_records)
    print("Number of gaps per sequence:")
    print(gaps_list)
    total_records -= 1

    new_gaps = 0
    new_count = 0
    x = 0
    conditional = True 

    for x in range(0, len(gaps_list)):
        print("Iterating over sequences|\tSequence: {}\tNumGaps:{} ".format( x,gaps_list[x]))
        if len(record_list[x].seq) >= num_percentile:    
            seq_lens.append(len(record_list[x].seq))
            seq_input = record_list[x].seq
            sequenceLength(seq_input) 
            GC_Content(seq_input) 
            beg_seq, end_seq, updated_seq = snip_seq(seq_input, gap_size) 
            poss_gaps=gaps_list[x] 
            if (len(updated_seq) > seq_size): #(num_gen <= probability)
                #poss_gaps = getNumGaps(num_gaps, ng_fixed)
                start_rand = getSequenceSize(seq_size, seq_fixed)
            
                gap_coordinates = validCoordinate(updated_seq, poss_gaps, start_rand)
                new_seq, seq_track = create_gaps(updated_seq, gap_coordinates, gap_size)

                all_gaps.append(seq_track)

                beginning  = "".join(beg_seq)
                end = "".join(end_seq)
                new = "".join(new_seq)

                final_seq = Seq(beginning+new+end)
                record = SeqRecord(seq=final_seq, id=record_list[x].id, name=record_list[x].name, description= record_list[x].description)

                if counts == 0:
                    mode = "w"
                else:
                    mode = "a"

                with open(args.output_file, mode) as output_handle:
                    count = SeqIO.write([record], output_handle, "fasta")

                for i in range(0, len(gap_coordinates)): 
                    gap_coordinates[i][0] += gap_size
                    gap_coordinates[i][1] += gap_size

                
                gap_coordinates.sort()

                coor_dict[record_list[x].id] = gap_coordinates
    
                fcount = flank_fasta(flank, flanks_fasta, seq_input, identifier, gap_coordinates, mode, fcount)
            

                    #append_gap_fasta(gap_fasta, coor_dict, seq_track)
                    # Generate a list of all the places that there aren't coordinates
                counts += 1
            else:
                record = SeqRecord(seq=seq_input, id=record_list[x].id, name=record_list[x].name, description= record_list[x].description)
                gap_coordinates = "None"
                coor_dict[record_list[x].id] = gap_coordinates

                if counts == 0:
                    mode = "w"
                else: 
                    mode = "a"
                with open(args.output_file, mode) as output_handle:
                    count = SeqIO.write([record], output_handle, "fasta")
                counts += 1
        x += 1

#        else:
#            continue
    
    ### Creation of the AGP file, gap_fasta file, and the bed file. 
    create_agp(file_input, agp_file, species)
    append_agp(coor_dict, agp_file, seq_lens, gap_size, species)
    create_gap_fasta(gap_fasta, coor_dict, all_gaps, identifier)
    create_bed(bed_file, coor_dict, identifier, flank)
    """
    This next portion of the program is where we will compute CPU usage, storage and memory usage, and other uses of this computers operating system 
    """

    cpu = os.times()
    ram = psutil.virtual_memory()
    disk = psutil.disk_usage(file_input)

    print(cpu)
    print()
    print(ram)
    print()
    print(disk)
    print()
