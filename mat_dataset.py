# ----------------------------------------------------------------------------------
# Functions to handle the unique attributes of the consensus1.mat file in Gong Lab
# @author: davidebear
# 
# Note: Refer to https://github.com/emiliosalazar/matToPython/blob/main/MatFileMethods.py for a more encompassing MatLabLoader
# ----------------------------------------------------------------------------------

# Needed imports
import numpy as np
import pandas as pd
import h5py

def load_matlab_file(mat_file, variable_name): # Note the double return statement
    """_summary_
    //UNIQUE: Only one variable for this .mat file
    
    Args:
        mat_file (string): Name of .mat file to intake (path must be given if not in current dir)
        variable_name (string): which variable from the mat file. Only one here because of the .mat file structure.

    Returns:
        h5py_object, data: Overall h5py_object and the specific variable as an nparray
    """
    h5py_object = h5py.File(mat_file, 'r')
    data = h5py_object.get(variable_name)
    data = np.array(data)
    return h5py_object, data

class DNA_SeqBlocks():
    """_summary_
    Loads the mat file with DNA sequences for easy iteration through all stored DNA sequence lists of varying coverage
    """
    def __init__(self, h5py_object, data):
        """_summary_
        Args:
            h5py_object (h5py_object): first return variable of load_matlab_file
            data (np.array): second return variable of load_matlab_file
        """
        self.h5py_object = h5py_object
        self.data = np.transpose(data)
        self.size = len(data) 
    def get_seqblock(self, number): # //TODO: make sure the default args work
        x = self.h5py_object[self.data[number-1, 0]][:, :] # //UNIQUE: 0th column because .mat file, align3 variable only has one row
        return np.transpose(x)
    
# CURRENT ISSUE:  I'm getting an int array instead of a char array. Use pd dataframe? How do I force a char?
def seqblock_parser(seqblock):
    """_summary_

    Args:
        seqblock (DNA_SeqBlocks.get_seqblock): A seqblock created from the DNA_SeqBlocks class

    Returns:
        CleanSeqBlock object : Divided aspects of a SeqBlock all in proper char format
    """
    seqblock_parsed = CleanSeqBlock()
    total_rows = len(np.transpose(seqblock)) # is there a more efficient way?
    total_columns = len(seqblock)
    number_of_reads = (total_rows - 4)/2 # last three rows have intreptations, first row is non-mutated target sequence, and N quality scores for N reads
    
    seqblock_parsed.reads_count = number_of_reads;
    seqblock_parsed.seq_len = total_rows;
    for i in range(total_rows):
        current_seq = seqblock[i]
        new_format = np.zeros((1,total_columns))
        for j in range(len(seqblock)):
            new_format[j] = chr(current_seq[j])
        if (i == 0):
            seqblock_parsed.target = new_format
        elif (i == total_rows - 1): # last row is mutations 'x', total_rows is one more than total index
            seqblock_parsed.interp_mutations = new_format
        elif (i == total_rows - 2): # 2nd last is subjective consensus
            seqblock_parsed.interp_consensus = new_format
        elif (i == total_rows - 3): 
            seqblock_parsed.interp_changes = new_format
        elif (i > 0 & i < number_of_reads + 1):
            seqblock_parsed.reads.append(new_format)
        else:
            seqblock_parsed.quality.append(new_format)
    return seqblock_parsed
            
    
class CleanSeqBlock():
    def __init__(self):
        self.reads_count = 0;
        self.seq_len = 0;
        
        self.target = 0;
        
        self.reads = [];
        self.quality = [];
        
        self.interp_changes = 0; # what is this again? 
        self.interp_consensus = 0;
        self.interp_mutations = 0;
    
