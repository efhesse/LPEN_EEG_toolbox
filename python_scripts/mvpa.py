import numpy as np
import scipy.io
import mne
from sklearn.metrics import roc_auc_score
from sklearn.cross_validation import StratifiedKFold
from mne.decoding import TimeDecoding
from mne.decoding import GeneralizationAcrossTime
import scipy.stats as spstat
import matplotlib.pyplot as plt

def create_mne_epochs(channels, data, events, sfreq,t_min,t_max,base):

    info = mne.create_info(channels, sfreq, ch_types='eeg', montage=None) 
                        
    if base:    
        epochs = mne.EpochsArray(data, info, events, tmin=t_min,tmax=t_max, event_id=None, reject=None, flat=None,
                        reject_tmin=None, reject_tmax=None, baseline=base, verbose=None)
    else:                        
        epochs = mne.EpochsArray(data, info, events, tmin=t_min,tmax=t_max, event_id=None, reject=None, flat=None,
                        reject_tmin=None, reject_tmax=None, baseline=None, verbose=None)                        
    
    return epochs
    
def create_y_data_for_decoding(epochs,event_1,event_0):
    events = epochs.events.astype(int)  
    #event_indexes = np.where((events[:,2] == event_1) | (events[:,2] == event_0))    
    #events = np.squeeze(events[event_indexes,:]) 
    allevents = event_1 + event_0
    events = events[np.logical_or.reduce([events[:,2] == x for x in allevents])]
    
    #get only relevant data trials
    data = epochs.get_data()    
    #data = np.squeeze(data[event_indexes,:,:])
    data = data[np.logical_or.reduce([events[:,2] == x for x in allevents])]
        
    y = np.zeros(len(events), dtype=int)
    #y[events[:,2] == event_1] = 1
    #y[events[:,2] == event_0] = 0
    y[np.where([events[:,2] == x for x in event_1])[1]] = 1
   
    info = mne.create_info(epochs.ch_names, epochs.info.get('sfreq'), ch_types='eeg', montage=None) 
    epochs_decode = mne.EpochsArray(data, info, events, tmin=epochs.tmin, event_id=None, reject=None, flat=None,
                        reject_tmin=None, reject_tmax=None, baseline=None, verbose=None)    

    return y,epochs_decode    
    
def calculate_diagonal_decoding(epochs,picked_channels,y,fold_nr,outputfile):    
    
    if picked_channels:
        epochs.pick_channels(picked_channels)
    else:
        print 'no channels to pick'

    #epochs.equalize_event_counts(epochs.event_id, copy=False)
    
    cv = StratifiedKFold(y=y,n_folds=fold_nr) # do a stratified cross-validation
    
    #armar clf y pasarlo al time decoding para mas de dos condiciones
    td = TimeDecoding(cv=cv, scorer=roc_auc_score, n_jobs=1,score_mode='fold-wise')
    
    # Fit, score, and plot
    td.fit(epochs, y=y)
    scores = td.score(epochs)
    
    if outputfile:
        scipy.io.savemat(outputfile,{'scores':scores})    
    print('TD DONE')
    return scores
    
def calculate_generalization_across_time_decoding(epochs,picked_channels,y,fold_nr,resamplefreq,outputfile):
   
        
    if picked_channels:
        epochs.pick_channels(picked_channels)
    else:
        print 'no channels to pick'
        
    if resamplefreq:
        #resample        
        epochs.resample(resamplefreq, npad=100, window='boxcar', n_jobs=1, copy=False, verbose=None)
        
    cv = StratifiedKFold(y=y,n_folds=fold_nr) # do a stratified cross-validation
    gat = GeneralizationAcrossTime(predict_mode='cross-validation', n_jobs=1, cv=cv, scorer=roc_auc_score,score_mode='fold-wise')

    # fit and score
    gat.fit(epochs, y=y)
    scores = gat.score(epochs)
    if outputfile:
        scipy.io.savemat(outputfile,{'scores':scores})    
        
    print('GAT DONE')
    return scores    
    
def plot_diagonal_decoding(epochs,scores,semfactor,outputfile):

    if type(scores) == list:
        scores = np.array(scores)
               
    times = epochs.times    
    chance_level = 0.5
    signlim = 0.3
    
    ci95 = spstat.sem(scores, axis=1) * semfactor
    scores_mean = np.mean(scores, axis=1)
    
    ci_lower = scores_mean - ci95
    ci_upper = scores_mean + ci95
       
    timepoints = scores_mean.shape[0]
        
    fig = plt.figure(figsize=(10, 10),dpi=300,edgecolor='white',facecolor='white')
    ax = plt.axes(frameon=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    chancelevel = [chance_level] * timepoints
        
    hscores, = plt.plot(times, scores_mean)
    hchance, = plt.plot(times,chancelevel,linestyle="--",color='k')
    
    plt.fill_between(times, ci_upper, ci_lower, color='dodgerblue', alpha='0.1')    
    
    #plot de indices significativos        
    statsig = ci_lower    
    statsig[statsig<= chance_level] = 0
    statsig[statsig!= 0] = 1
    statsig_final = statsig[statsig==1]*signlim
    tstats = times[statsig==1]
    
    plt.plot(tstats, statsig_final, marker='.',markersize=0.7,linestyle='None')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Classification score (AUC)')       
      
    plt.subplots_adjust(top=0.78)
    plt.show()    
    fig.tight_layout()  
    
    if outputfile:
        plt.savefig(outputfile + '.png')    
        #plt.savefig(outputfile + '.eps' , format='eps', dpi=300)    
    
    print('DONE saving TD plot')    
    
def plot_generalization_across_time(epochs,scores,semfactor,outputfile):
    
    chance_level = 0.5    

    ci95 = spstat.sem(scores, axis=2) * semfactor
    scores_mean = np.mean(scores, axis=2)
    
    ci_lower = scores_mean - ci95            
    times = epochs.times    

    # Matriz
    #fig = plt.figure()        
    fig, ax = plt.subplots(1, 1, figsize=(16, 16),dpi=300,edgecolor='white',facecolor='white')
    sig_mask = ci_lower > chance_level
    plt.imshow(
        np.ma.masked_array(scores_mean, mask=sig_mask),
        interpolation='nearest', origin='lower',
        extent=[times[0], times[-1], times[0], times[-1]],
        cmap='gray')
    vmin=0.5
    vmax=1
    plt.imshow(
        np.ma.masked_array(scores_mean, mask=~sig_mask),
        interpolation='nearest', origin='lower',
        extent=[times[0], times[-1], times[0], times[-1]],vmin=vmin, vmax=vmax,
        cmap='RdBu_r')

    if min(times) <= 0:
        plt.axvline(0, color='k')
        plt.axhline(0, color='k')
        
    #plt.colorbar(ticks = colorbarticks)
    plt.colorbar()
    
    plt.xlabel('Testing time (ms)')
    plt.ylabel('Training time (ms)')   

    ax.axvline(0, color='k')
    ax.axhline(0, color='k')
    
    plt.show()
    fig.tight_layout()    
        
    if outputfile:
        plt.savefig(outputfile + '.png')    
        #plt.savefig(outputfile + '.eps' , format='eps', dpi=300)    
      
    print('DONE saving GAT plot')
    