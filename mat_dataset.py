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

class DNA_Sequences_Dataset():
    def __init__(self, h5py_object, data):
        self.h5py_object = h5py_object
        self.data = data
        self.size = len(data.shape)
    def get_sequences(self, number): # //TODO: make sure the default args work
        x = self.h5py_object[self.data[number-1, 0]][:, :] # //UNIQUE: 0th column because .mat file, align3 variable only has one row
        return x
    
# CURRENT ISSUE:  I'm getting an int array instead of a char array. Use pd dataframe? How do I force a char?

