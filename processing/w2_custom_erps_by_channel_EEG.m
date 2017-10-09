function [data] = w2_custom_erps_by_channel_EEG(path_to_files,condition_1,condition_2,cycles,freq_range,alpha,fdr,scale,basenorm,erps_max,path_to_save,data);

%-------DESCRIPTION---------------
%SEE NEW_TIMEF OR w2_erps_by_channel
%----------------------------------

%-------PATH MANAGEMENT-----------
%----------------------------------
%add path of preprocessing 
addpath(genpath(fullfile(data.ieeglab_path,'processing')));

%modifies paths to include parent directory
path_to_files = fullfile(data.parent_directory, path_to_files);
path_to_save = fullfile(data.parent_directory, path_to_save);

%check if path to files exist
assert(exist(path_to_files, 'dir') == 7,['Directory not found! path_to_files: ' path_to_files '.']);

%create directory where the trimmed sets will be stored
if ~exist(path_to_save, 'dir')
  mkdir(path_to_save);
end

%-------LOAD PARAMETERS------------
%----------------------------------
cycles = str2num(cycles);
freq_range = str2num(freq_range);
alpha = str2num(alpha);
if alpha == 0
    alpha = NaN;
end
erps_max = str2num(erps_max);

%---------RUN---------------------
%---------------------------------

%load .set 
files = dir(fullfile(path_to_files,'*.set'));
filenames = {files.name}';  
file_nr = size(filenames,1);

for suj = 1 : file_nr
    file_name = filenames{suj};
    disp(path_to_files)
    disp(file_name)
    
    channel_labels = [];
    %4D matrix where results will be stored (channels, freq, times, subjects)    
    [c1_erps,c2_erps,c1_c2_erps] = deal([],[],[]);
    %5D matrix single trial results (channels,freq,times,epochs,subjects)
    [c1_tfX,c2_tfX] = deal([],[]);
    %4D matrix where statistical results will be stored (channels, freq,
    %times, subjects)
    [c1_erpsboot,c2_erpsboot,c1_c2_erpsboot] = deal([],[],[]);
    %2D matrix, where baseline powers are stored (base,subjects)
    [c1_mbases,c2_mbases,c1_c2_mbases] = deal([],[],[]);
    [c1_data,c2_data] = deal([],[]);
    [c1_itc,c2_itc,c1_c2_itc] = deal([],[],[]);
    [c1_itcboot,c2_itcboot,c1_c2_itcboot] = deal([],[],[]);
    [c1_maskersp,c2_maskersp] = deal([],[]);
    [c1_maskitc,c2_maskitc] = deal([],[]);
    [resdiff1,resdiff2] = deal([],[]);
    [c1_pa,c2_pa] = deal([],[]);
    
    %load set
    EEG = pop_loadset('filename',file_name,'filepath', path_to_files);
    ch_nr = length(EEG.chanlocs); %number of channel
    
    %load EEG according to condition selections
    two_conditions = 0;
    c2_EEG = [];
    if ~isempty(condition_1) && ~isempty(condition_2) %two conditions
        two_conditions = 1;
        c1_EEG = pop_selectevent( EEG, 'type', {condition_1} ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        c2_EEG = pop_selectevent( EEG, 'type', {condition_2} ,'deleteevents','off','deleteepochs','on','invertepochs','off');
    elseif isempty(condition_1) && isempty(condition_2) %no filtering by condition, one dataset
        c1_EEG = EEG;
    else    %only one condition filtering, one dataset
        if isempty(condition_1) && ~isempty(condition_2) %
            condition_1 = condition_2;
            condition_2 = '';
            c1_EEG = pop_selectevent( EEG, 'type', {condition_2} ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        else
            c1_EEG = pop_selectevent( EEG, 'type', {condition_1} ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        end
    end 
        
    %calculate point range 
    tlimits = [EEG.xmin, EEG.xmax]*1000;
    pointrange1 = round(max((tlimits(1)/1000-EEG.xmin)*EEG.srate, 1));
    pointrange2 = round(min((tlimits(2)/1000-EEG.xmin)*EEG.srate, EEG.pnts));
    pointrange = [pointrange1:pointrange2];
     
    for ch = 1 : 3        
    %for ch = 1 : ch_nr        
        chanlabel = EEG.chanlocs(ch).labels;   
        channel_labels{ch} = chanlabel;
        
        c1_tmpsig = c1_EEG.data(ch,pointrange,:);
        c1_tmpsig = reshape( c1_tmpsig, length(ch), size(c1_tmpsig,2)*size(c1_tmpsig,3));
        data_to_process = c1_tmpsig(:,:);
        
        if ~isempty(condition_2)
            c2_tmpsig = c2_EEG.data(ch,pointrange,:);
            c2_tmpsig = reshape( c2_tmpsig, length(ch), size(c2_tmpsig,2)*size(c2_tmpsig,3));
            data_to_process = {c1_tmpsig(:,:),c2_tmpsig(:,:)};
        end

        %calculate timefreq 
        %[P,R,mbase,timesout,freqs,Pboot,Rboot,alltfX,g] = newtimef_2_conditions( {c1_tmpsig(:, :),c2_tmpsig(:, :)}, length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles, 'plotersp','off', 'plotitc' , 'off','topovec', ch, 'elocs', EEG.chanlocs,'title',{condition_1 condition_2},'freqs',freq_range,'alpha',alpha,'mcorrect',fdr,'scale',scale,'basenorm',basenorm,'erspmax',erps_max, 'ntimesout', 400, 'padratio', 4,'baseline',[0],'caption',chanlabel) ;                        
        if two_conditions
            %[P,R,mbase,timesout,freqs,Pboot,Rboot,alltfX,g] = newtimef_2_conditions( data_to_process, length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles, 'plotersp','off', 'plotitc' , 'off','topovec', ch, 'elocs', EEG.chanlocs,'title',{condition_1 condition_2},'freqs',freq_range,'alpha',alpha,'mcorrect',fdr,'scale',scale,'basenorm',basenorm,'erspmax',erps_max, 'ntimesout', 400, 'padratio', 4,'baseline',[0],'caption',chanlabel) ; 
            [ERP,P,R,mbase,timesout,freqs,Pboot,Rboot,resdiff,alltfX,PA,maskersp,maskitc,g] = custom_newtimef(data_to_process, length(pointrange), tlimits, EEG.srate, cycles,[], 'plotersp','off', 'plotitc' , 'off','topovec', ch, 'elocs', EEG.chanlocs,'title',{condition_1 condition_2},'freqs',freq_range,'alpha',alpha,'mcorrect',fdr,'scale',scale,'basenorm',basenorm,'erspmax',erps_max, 'ntimesout', 400, 'padratio', 4,'baseline',[0],'caption',chanlabel);
        else
            %[P,R,mbase,timesout,freqs,Pboot,Rboot,alltfX,PA,g] = m_newtimef( data_to_process, length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles, 'plotersp','off', 'plotitc' , 'off','topovec', ch, 'elocs', EEG.chanlocs,'title',{condition_1 condition_2},'freqs',freq_range,'alpha',alpha,'mcorrect',fdr,'scale',scale,'basenorm',basenorm,'erspmax',erps_max, 'ntimesout', 400, 'padratio', 4,'baseline',[0],'caption',chanlabel) ;                        
            [ERP,P,R,mbase,timesout,freqs,Pboot,Rboot,resdiff,alltfX,PA,maskersp,maskitc,g] = custom_newtimef(data_to_process, length(pointrange), tlimits, EEG.srate, cycles,[], 'plotersp','off', 'plotitc' , 'off','topovec', ch, 'elocs', EEG.chanlocs,'title',condition_1,'freqs',freq_range,'alpha',alpha,'mcorrect',fdr,'scale',scale,'basenorm',basenorm,'erspmax',erps_max, 'ntimesout', 400, 'padratio', 4,'baseline',[0],'caption',chanlabel);
        end        
        
        %load results in matrices of overall results            
        if two_conditions
            c1_data(ch,:,suj) = ERP{1};
            c2_data(ch,:,suj) = ERP{2};
            c1_itc(ch,:,:,suj) = R{1}; 
            c2_itc(ch,:,:,suj) = R{2}; 
            c1_c2_itc(ch,:,:,suj) = R{3}; 
            c1_erps(ch,:,:,suj) = P{1};
            c2_erps(ch,:,:,suj) = P{2};
            c1_c2_erps(ch,:,:,suj) = P{3};

            if ~isnan(alpha)
                c1_erpsboot(ch,:,:,suj) = Pboot{1};
                c2_erpsboot(ch,:,:,suj) = Pboot{2};
                c1_c2_erpsboot(ch,:,:,:,suj) = Pboot{3};
                c1_itcboot(ch,:,suj) = Rboot{1};
                c2_itcboot(ch,:,suj) = Rboot{2};
                c1_c2_itcboot(ch,:,:,suj) = Rboot{3};
                c1_maskersp(ch,:,:,suj) = maskersp{1};
                c2_maskersp(ch,:,:,suj) = maskersp{2};
                c1_maskitc(ch,:,:,suj) = maskitc{1};
                c2_maskitc(ch,:,:,suj) = maskitc{2};
            end

            c1_tfX(ch,:,:,:) = alltfX{1};
            c2_tfX(ch,:,:,:) = alltfX{2};
            c1_mbases(ch,:,suj) = mbase{1};        
            c2_mbases(ch,:,suj) = mbase{2}; 
            c1_c2_mbases(ch,:,suj) = mbase{3};
            c1_pa(ch,:,:,:) = PA{1}; 
            c2_pa(ch,:,:,:) = PA{2}; 
            resdiff1(ch,:,:,suj) = resdiff{1}; 
            resdiff2(ch,:,:,suj) = resdiff{2};
        else
            c1_data(ch,:,suj) = ERP;
            c1_itc(ch,:,:,suj) = R; 
            c1_erps(ch,:,:,suj) = P;
            if ~isnan(alpha)
                c1_erpsboot(ch,:,:,suj) = Pboot;   
                c1_itcboot(ch,:,suj) = Rboot;
                c1_maskersp(ch,:,:,suj) = maskersp;
                c1_maskitc(ch,:,:,suj) = maskitc;
            end
            c1_tfX(ch,:,:,:) = alltfX;
            c1_mbases(ch,:,suj) = mbase;        
            c1_pa(ch,:,:,:) = PA; 
        end
    end
    
    %save data for suj
    [filepath,file_name_to_save,ext] = fileparts(file_name);
    if two_conditions
        mat_name = fullfile(path_to_save,[file_name_to_save '_' condition_1 '_' condition_2 '.mat']);
        s_data = {c1_data(:,:,suj),c2_data(:,:,suj)};
        s_erps = {c1_erps(:,:,:,suj),c2_erps(:,:,:,suj),c1_c2_erps(:,:,:,suj)};
        s_erpsboot = {c1_erpsboot(:,:,:,suj),c2_erpsboot(:,:,:,suj),c1_c2_erpsboot(:,:,:,suj)};
        s_itc = {c1_itc(:,:,suj),c2_itc(:,:,suj),c1_c2_itc(:,:,suj)};
        s_itcboot = {c1_itcboot(:,:,suj),c2_itcboot(:,:,suj),c1_c2_itcboot(:,:,suj)};
        s_tfX = {c1_tfX,c2_tfX};
        s_resdiff = {resdiff1(:,:,suj),resdiff2(:,:,suj)};
        s_mbases = {c1_mbases(:,:,suj),c2_mbases(:,:,suj),c1_c2_mbases(:,:,suj)};
        s_maskerps = {c1_maskersp(:,:,suj),c2_maskersp(:,:,suj)};
        s_maskitc = {c1_maskitc(:,:,suj), c2_maskitc(:,:,suj)};
        s_pa = {c1_pa,c2_pa};
        
        c1_alltfX(suj).tfX = c1_tfX;
        c2_alltfX(suj).tfX = c2_tfX;
        c1_allpa(suj).PA =  c1_pa;
        c2_allpa(suj).PA = c2_pa;
    else
        mat_name = fullfile(path_to_save,[file_name_to_save '_' condition_1 '.mat']);
        s_resdiff = [];
        s_data = c1_data(:,:,suj);
        s_itc = c1_itc(:,:,suj);
        s_erps = c1_erps(:,:,:,suj);
        s_erpsboot = c1_erpsboot(:,:,:,suj);
        s_itcboot = c1_itcboot(:,:,suj);
        s_maskerps = c1_maskersp(:,:,suj);
        s_maskitc = c1_maskitc(:,:,suj);
        s_tfX = c1_tfX;
        s_mbases = c1_mbases(:,:,suj);    
        s_pa = c1_pa;
        
        c1_alltfX(suj).tfX = c1_tfX;
        c1_allpa(suj).PA =  c1_pa; 
    end

    save(mat_name, 's_erps','s_erpsboot','s_tfX','s_mbases','s_resdiff','s_data','s_itc','s_itcboot','s_maskerps','s_maskitc','s_pa','timesout','freqs','g','channel_labels');
    
    clear s_erps s_erpsboot s_tfX s_mbases s_resdiff s_data s_itc s_itcboot s_maskerps s_maskitc s_pa
end
%save results
if two_conditions
    erps = {c1_erps,c2_erps,c1_c2_erps};
    erpsboot = {c1_erpsboot,c2_erpsboot,c1_c2_erpsboot};
    itc = {c1_itc, c2_itc, c1_c2_itc}; 
    itcboot = {c1_itcboot,c2_itcboot,c1_c2_itcboot};
    resdiff = {resdiff1,resdiff2};
    mdata = {c1_data,c2_data};
    tfX = {c1_alltfX,c2_alltfX};
    pa = {c1_allpa,c2_allpa};
    mbases = {c1_mbases,c2_mbases,c1_c2_mbases};
    maskerps = {c1_maskersp,c2_maskersp};
    maskitc = {c1_maskitc,c2_maskitc};
    prefix_file_name_to_save = [condition_1 '_' condition_2];
else
    erps = c1_erps;
    erpsboot = c1_erpsboot;
    mdata = c1_data;
    itc = c1_itc;
    itcboot = c1_itcboot;
    tfX = c1_alltfX;
    pa = c1_allpa;
    mbases = c1_mbases;
    resdiff = [];
    maskerps = c1_maskersp;
    maskitc = c1_maskitc;
    prefix_file_name_to_save = [condition_1];
end

%save results in mat
mat_name = fullfile(path_to_save,[prefix_file_name_to_save '.mat']);
save(mat_name, 'erps','erpsboot','tfX','mbases','timesout','freqs','g','mdata','itc','itcboot','resdiff','maskerps','maskitc','pa','channel_labels');
