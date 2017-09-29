# -*- coding: utf-8 -*-
"""
Created on Mon Oct 03 17:26:02 2016

@author: Ineco
"""

import sys
import mne
import scipy.io
from mvpa import *

matlab_args = sys.argv[2]
input_filename = sys.argv[1]#ESTE SET TIENE QUE TENER LAS CONDICIONES CODIFICADAS POR NUMEROS -
    #VUELVO A LA IDEA DE QUIZAS DENTRO DE TU ESTRUCTURA TENER UN DICCIONARIO PARA QUE EL USUARIO PUEDA PONER INT Y VOS INTERPRETARLO COMO 1
    #PORQUE SIEMPRE ES MAS FACIL DE ENTENDER ASI --- NO SE QUE OPINAS?
mat = scipy.io.loadmat(matlab_args, chars_as_strings=True)
python_args = mat['python_args']

#PARAMETERS
def get_array(array):
    if len(array[0][0]) > 0:
        return array[0][0][0]
    return []

picked_channels = get_array(python_args['picked_channels'])
event_1 = get_array(python_args['event_1'])
event_0 = get_array(python_args['event_0'])

fold_nr = python_args['fold_nr'][0][0][0][0]
outputfile = python_args['outputfile'][0][0][0]
semfactor = python_args['semfactor'][0][0][0][0]

epochs = mne.read_epochs_eeglab(input_filename)
y,epochs_decode = create_y_data_for_decoding(epochs,event_1,event_0)

scores = calculate_diagonal_decoding(epochs_decode,picked_channels,y,fold_nr,outputfile)
plot_diagonal_decoding(epochs_decode,scores,semfactor,outputfile)