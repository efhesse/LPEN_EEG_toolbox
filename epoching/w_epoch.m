function [ EEG, data ] = w_epoch(epoch_window, base_window,file_to_save, EEG, data)
types = unique({EEG.event.type});
str = '';
for i=1:length(types)
    if isequal(str,'')
        str = types{i};
    else
        str = [str ' ' types{i}];
    end
end
answer = inputdlg('Types','Types', 1, {str});
if ~isempty(answer)
    types = strsplit(answer{1}, ' ');
end
epoch_window = str2num(epoch_window);
base_window = str2num(base_window);

if ~isequal(file_to_save,'')
    file_to_save = fullfile(data.path, [file_to_save '.mat']);
end

EEG = ieeg_epoch(types, epoch_window, base_window,file_to_save, EEG);

end

