from mat_dataset import *

file_name = 'consensus1.mat'
variable = 'align3'

h5py_file, data = load_matlab_file(file_name, variable)

Sequence_1009 = DNA_SeqBlocks(h5py_file, data);