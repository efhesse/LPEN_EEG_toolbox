import scipy.io
import mne
import sys
from mvpa import *
args = sys.argv[2]
input_fname = sys.argv[1]
print args
mat = scipy.io.loadmat(args, chars_as_strings=True)
print mat
mat = mat['python_args']
#one_value = mat['one_value'][0][0][0][0]
#double_array = mat['double_array'][0][0][0]
#string_value = mat['path'][0][0][0]
#print one_value, double_array, string_value

def get_array(array):
    if len(array[0][0]) > 0:
        return array[0][0][0]
    return []
picked_channels = get_array(mat['picked_channels'])
event_1 = get_array(mat['event_1'])
event_0 = get_array(mat['event_0'])

fold_nr = mat['fold_nr'][0][0][0][0]
outputfile = mat['outputfile'][0][0][0]
semfactor = mat['semfactor'][0][0][0][0]

#import ipdb; ipdb.set_trace()
epochs = mne.read_epochs_eeglab(input_fname)
y,epochs_decode = create_y_data_for_decoding(epochs,event_1,event_0)

scores = calculate_diagonal_decoding(epochs_decode,picked_channels,y,fold_nr,outputfile)
plot_diagonal_decoding(epochs_decode,scores,semfactor,outputfile)
#eeg = mat['eeg_data']
#eeg = eeg*3
#scipy.io.savemat(matrix, {'eeg_data': eeg})