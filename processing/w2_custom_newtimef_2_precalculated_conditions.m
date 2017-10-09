function w2_custom_newtimef_2_precalculated_conditions(path_to_files,condition1_data,condition2_data,path_to_save)

%[ERP,P,R,mbase,timesout,freqs,Pboot,Rboot,resdiff,alltfX,PA,maskersp,maskitc,g] =
clear freqs
%load condition 1
c1 = load(fullfile(path_to_files,condition1_data));
all_ERP1 = c1.mdata;
all_P1 = c1.erps;
all_R1 = c1.itc;
all_mbase1 = c1.mbases;
all_Pboot1 = c1.erpsboot;
all_Rboot1 = c1.itcboot;
channel_nr = length(c1.channel_labels);
channel_labels = c1.channel_labels;
assert(channel_nr == size(all_P1,1),'Error. Number of channel labels do not match calculated channel number.');
suj_nr = size(all_P1,4);
%resdiff -> should be empty, will be calculated in this function
all_tfX1 = c1.tfX; %aca seguro hay que hacer algo

all_PA1 = c1.pa; 
all_maskersp1 = c1.maskerps;
all_maskitc1 = c1.maskitc;
g1 = c1.g;
condition_1 = g1.title;

%load condition 2
c2 = load(fullfile(path_to_files,condition2_data));
all_ERP2 = c2.mdata;
all_P2 = c2.erps;

assert(channel_nr == size(all_P2,1),'Error. Different channel number for conditions.');
assert(suj_nr == size(all_P2,4),'Error. Different subject number for conditions.');

all_R2 = c2.itc;
all_mbase2 = c2.mbases;
all_Pboot2 = c2.erpsboot;
all_Rboot2 = c2.itcboot;
%resdiff -> should be empty, will be calculated in this function
all_tfX2 = c2.tfX;

all_PA2 = c2.pa; 
all_maskersp2 = c2.maskerps;
all_maskitc2 = c2.maskitc;
g2 = c2.g;
condition_2 = g2.title;

%result matrices
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

%mean over sujbects by channel 
for ch = 1 : channel_nr
    suj.ERP1 = squeeze(mean(all_ERP1(ch,:,:),3));
    suj.ERP2 = squeeze(mean(all_ERP2(ch,:,:),3));
    suj.P1 = squeeze(mean(all_P1(ch,:,:,:),4));
    suj.P2 = squeeze(mean(all_P2(ch,:,:,:),4));
    suj.R1 = squeeze(mean(all_R1(ch,:,:,:),4));
    suj.R2 = squeeze(mean(all_R2(ch,:,:,:),4));
    suj.mbase1 = squeeze(mean(all_mbase1(ch,:,:),3));
    suj.mbase2 = squeeze(mean(all_mbase2(ch,:,:),3));
    suj.Pboot1 = squeeze(mean(all_Pboot1(ch,:,:,:),4));
    suj.Pboot2 = squeeze(mean(all_Pboot2(ch,:,:,:),4));
    suj.Rboot1 = squeeze(mean(all_Rboot1(ch,:,:),3));
    suj.Rboot2 = squeeze(mean(all_Rboot2(ch,:,:),3));

    for i = 1 : suj_nr
        temp_alltfX1(:,:,i) = squeeze(mean(all_tfX1(i).tfX(ch,:,:,:),4));
        temp_alltfX2(:,:,i) = squeeze(mean(all_tfX2(i).tfX(ch,:,:,:),4));
        temp_allpa1(:,:,i) = squeeze(mean(all_PA1(i).PA(ch,:,:,:),4));
        temp_allpa2(:,:,i) = squeeze(mean(all_PA2(i).PA(ch,:,:,:),4));
    end

    suj.alltfX1 = temp_alltfX1;
    suj.alltfX2 = temp_alltfX2;
    suj.PA1 = temp_allpa1;
    suj.PA2 = temp_allpa2;
    suj.maskersp1 = squeeze(mean(all_maskersp1(ch,:,:,:),4));
    suj.maskersp2 = squeeze(mean(all_maskersp2(ch,:,:,:),4));
    suj.maskitc1 = squeeze(mean(all_maskitc1(ch,:,:,:),4));
    suj.maskitc2 = squeeze(mean(all_maskitc2(ch,:,:,:),4));
    
    suj.timesout = c1.timesout;
    suj.freqs = c1.freqs;
    c1.g.caption = c1.channel_labels(ch);
    c1.g.title = {condition_1,condition_2};
    suj.g = c1.g; %set other conditions such as alpha and fdr, for example, before calling precalculated function

    [ERP,P,R,mbase,timesout,freqsout,Pboot,Rboot,resdiff,alltfX,PA,maskersp,maskitc,g] = custom_newtimef_2_precalculated_conditions(suj);
    
    %load in results matrix
    c1_data(ch,:) = ERP{1};
    c2_data(ch,:) = ERP{2};
    c1_itc(ch,:,:) = R{1}; 
    c2_itc(ch,:,:) = R{2}; 
    c1_c2_itc(ch,:,:) = R{3}; 
    c1_erps(ch,:,:) = P{1};
    c2_erps(ch,:,:) = P{2};
    c1_c2_erps(ch,:,:) = P{3};

    if ~isnan(suj.g.alpha)
        c1_erpsboot(ch,:,:) = Pboot{1};
        c2_erpsboot(ch,:,:) = Pboot{2};
        c1_c2_erpsboot(ch,:,:,:) = Pboot{3};
        c1_itcboot(ch,:) = Rboot{1};
        c2_itcboot(ch,:) = Rboot{2};
        c1_c2_itcboot(ch,:,:,:) = Rboot{3};
        c1_maskersp(ch,:,:) = maskersp{1};
        c2_maskersp(ch,:,:) = maskersp{2};
        c1_maskitc(ch,:,:) = maskitc{1};
        c2_maskitc(ch,:,:) = maskitc{2};
    end

    c1_tfX(ch,:,:,:) = alltfX{1};
    c2_tfX(ch,:,:,:) = alltfX{2};
    c1_mbases(ch,:) = mbase{1};        
    c2_mbases(ch,:) = mbase{2}; 
    c1_c2_mbases(ch,:) = mbase{3};
    c1_pa(ch,:,:,:) = PA{1}; 
    c2_pa(ch,:,:,:) = PA{2}; 
    resdiff1(ch,:,:) = resdiff{1}; 
    resdiff2(ch,:,:) = resdiff{2};
end

%save results matrix
erps = {c1_erps,c2_erps,c1_c2_erps};
erpsboot = {c1_erpsboot,c2_erpsboot,c1_c2_erpsboot};
itc = {c1_itc, c2_itc, c1_c2_itc}; 
itcboot = {c1_itcboot,c2_itcboot,c1_c2_itcboot};
resdiff = {resdiff1,resdiff2};
mdata = {c1_data,c2_data};
tfX = {c1_tfX,c2_tfX};
pa = {c1_pa,c2_pa};
mbases = {c1_mbases,c2_mbases,c1_c2_mbases};
maskerps = {c1_maskersp,c2_maskersp};
maskitc = {c1_maskitc,c2_maskitc};
prefix_file_name_to_save = [condition_1 '_' condition_2];

mat_name = fullfile(path_to_save,[prefix_file_name_to_save '.mat']);
save(mat_name, 'erps','erpsboot','tfX','mbases','timesout','freqsout','g','mdata','itc','itcboot','resdiff','maskerps','maskitc','pa','channel_labels');