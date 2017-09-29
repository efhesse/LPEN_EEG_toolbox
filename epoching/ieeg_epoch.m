function [ EEG ] = ieeg_epoch(types, epoch_window, base_window,file_to_save, EEG)

    EEG = pop_epoch(EEG, types, epoch_window, 'newname', '', 'epochinfo', 'yes');

    if epoch_window(1) ~= 0
        display('Baseline Removal')
        EEG = pop_rmbase(EEG, base_window);
    end

    if ~isequal(file_to_save,'')
        data = EEG.data;
        eval(['save ' file_to_save ' data;'])
    end
end

