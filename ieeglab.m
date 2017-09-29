function varargout = ieeglab(varargin)
% IEEGLAB MATLAB code for ieeglab.fig
%      IEEGLAB, by itself, creates a new_project IEEGLAB or raises the existing
%      singleton*.
%
%      H = IEEGLAB returns the handle to a new_project IEEGLAB or the handle to
%      the existing singleton*.
%
%      IEEGLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IEEGLAB.M with the given input arguments.
%
%      IEEGLAB('Property','Value',...) creates a new_project IEEGLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ieeglab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ieeglab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ieeglab

% Last Modified by GUIDE v2.5 24-Oct-2016 12:25:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ieeglab_OpeningFcn, ...
                   'gui_OutputFcn',  @ieeglab_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ieeglab is made visible.
function ieeglab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ieeglab (see VARARGIN)
% addpath(genpath('/Users/maximo/Documents/MATLAB/ieeg/Scripts'));

data = struct();

data.data_status = 0;
[pathstr,~,~] = fileparts(which('ieeglab'));
data.ieeglab_path = pathstr;
%LOAD preprocessing
data.preprocessing_functions = cell(0);
[str, input] = read_file_by_lines(fullfile(data.ieeglab_path,'preprocessing','functions.txt'));
data.preprocessing_input = input;
set(handles.preprocessing_function_select_menu, 'String', str);
[data.current_preprocessing_function.str, data.current_preprocessing_function.pos] = get_current_popup_string(handles.preprocessing_function_select_menu);

%LOAD epoching
[str, input] = read_file_by_lines(fullfile(data.ieeglab_path,'epoching','functions.txt'));
set(handles.epoching_select_menu, 'String', str);
data.epoching_input = input;
[data.epoching_function.str, data.epoching_function.pos] = get_current_popup_string(handles.epoching_select_menu);

%LOAD processing
data.processing_functions = cell(0);
[str, input] = read_file_by_lines(fullfile(data.ieeglab_path,'processing','functions.txt'));
data.processing_input = input;
set(handles.processing_select_menu, 'String', str);
[data.current_processing_function.str, data.current_processing_function.pos] = get_current_popup_string(handles.processing_select_menu);

set(handles.figure1,'CloseRequestFcn',@close_GUI);

handles.data = data;
load_handles(handles);
% Choose default command line output for ieeglab
handles.output = hObject;
 
% Update handles structure
guidata(hObject, handles);
 
% UIWAIT makes ieeglab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ieeglab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function new_project_Callback(hObject, eventdata, handles)
[file_name, path] = uiputfile('*.mat','Save Workspace As');
if file_name
    handles.data.file_name = file_name;
    handles.data.path = path;
    save_file( handles);
    load_handles(handles);
    guidata(hObject,handles)
end

% --- Executes on selection change in preprocessing_function_select_menu.
function preprocessing_function_select_menu_Callback(hObject, eventdata, handles)
[handles.data.current_preprocessing_function.str, handles.data.current_preprocessing_function.pos] = get_current_popup_string(hObject);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function preprocessing_function_select_menu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject,handles)


% --- Executes on button press in add_button.
function add_button_Callback(hObject, eventdata, handles)
% hObject    handle to add_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
preprocessing_function.str = handles.data.current_preprocessing_function.str;
preprocessing_function.pos = handles.data.current_preprocessing_function.pos;
preprocessing_function.params = {};
input = handles.data.preprocessing_input{handles.data.current_preprocessing_function.pos};
if ~isempty(input(:,1))
    answer = dynamic_inputdlg(input(:,1),'Input', 1, input(:,2));
    if ~isempty(answer)
        preprocessing_function.params = answer;
    else
        return
    end
end
handles.data.preprocessing_functions{end+1} = preprocessing_function;
handles.data.data_status = 0;
load_handles(handles);
guidata(hObject,handles)

% --- Executes on selection change in pre_processing_functions_list.
function pre_processing_functions_list_Callback(hObject, eventdata, handles)
handles.data.current_preprocessing_function_to_delete = get(hObject,'Value');
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function pre_processing_functions_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in edit_preprocessing_btn.
function edit_preprocessing_btn_Callback(hObject, eventdata, handles)
if isfield(handles.data,'current_preprocessing_function_to_delete')
    preprocessing_function = handles.data.preprocessing_functions{handles.data.current_preprocessing_function_to_delete};
    input = handles.data.preprocessing_input{preprocessing_function.pos};
    if ~isempty(input(:,1))
        answer = dynamic_inputdlg(input(:,1),'Input', 1, preprocessing_function.params);
        if ~isempty(answer)
            preprocessing_function.params = answer;
            handles.data.preprocessing_functions{handles.data.current_preprocessing_function_to_delete} = preprocessing_function;
            handles.data.data_status = 0;
            load_handles(handles);
            guidata(hObject,handles)
        end
    end
end
% --- Executes on button press in delete_button.
function delete_button_Callback(hObject, eventdata, handles)
if isfield(handles.data,'current_preprocessing_function_to_delete')
    if ~isempty(handles.data.preprocessing_functions(:))
        handles.data.preprocessing_functions(:, handles.data.current_preprocessing_function_to_delete) = [];
        handles.data.data_status = 0;
    end
    load_handles(handles);
    guidata(hObject,handles)
end


% --- Executes on selection change in epoching_select_menu.
function epoching_select_menu_Callback(hObject, eventdata, handles)
[handles.data.epoching_function.str, handles.data.epoching_function.pos] = get_current_popup_string(hObject);
handles.data.data_status = min(1,handles.data.data_status);
handles.data.epoching_function.params = {};
input = handles.data.epoching_input{handles.data.epoching_function.pos};
handles.data.epoching_function.input = input;
if ~isempty(input(:,1))
    answer = dynamic_inputdlg(input(:,1),'Input', 1, input(:,2));
    if ~isempty(answer)
        handles.data.epoching_function.params = answer;
    else
        return
    end
end
load_handles(handles);
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function epoching_select_menu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% Save the handles structure.
guidata(hObject,handles)

% --- Executes on button press in edit_epoching_btn.
function edit_epoching_btn_Callback(hObject, eventdata, handles)
input = handles.data.epoching_input{handles.data.epoching_function.pos};
params = input(:,2);
if isfield(handles.data.epoching_function,'params')
    params = handles.data.epoching_function.params;
end
if ~isempty(input(:,1))
    answer = dynamic_inputdlg(input(:,1),'Input', 1, params);
    if ~isempty(answer)
        handles.data.epoching_function.params = answer;
        handles.data.data_status = min(1,handles.data.data_status);
        load_handles(handles);
        guidata(hObject,handles)
    end
end
% Save the handles structure.



% --- Executes on button press in edit_processing_btn.
function edit_processing_btn_Callback(hObject, eventdata, handles)
if isfield(handles.data,'current_processing_function_to_delete')
    processing_function = handles.data.processing_functions{handles.data.current_processing_function_to_delete};
    input = handles.data.processing_input{processing_function.pos};
    if ~isempty(input(:,1))
        answer = dynamic_inputdlg(input(:,1),'Input', 1, processing_function.params);
        if ~isempty(answer)
            processing_function.params = answer;
            handles.data.processing_functions{handles.data.current_processing_function_to_delete} = processing_function;
            handles.data.data_status = min(2,handles.data.data_status);
            load_handles(handles);
            guidata(hObject,handles)
        end
    end
end

% --- Executes on selection change in processing_select_menu.
function processing_select_menu_Callback(hObject, eventdata, handles)
[handles.data.current_processing_function.str, handles.data.current_processing_function.pos] = get_current_popup_string(hObject);
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function processing_select_menu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in add_processing_function_button.
function add_processing_function_button_Callback(hObject, eventdata, handles)
processing_function.str = handles.data.current_processing_function.str;
processing_function.pos = handles.data.current_processing_function.pos;
processing_function.params = {};
input = handles.data.processing_input{handles.data.current_processing_function.pos};
if ~isempty(input(:,1))
    answer = dynamic_inputdlg(input(:,1),'Input', 1, input(:,2));
    if ~isempty(answer)
        processing_function.params = answer;
    else
        return
    end
end
handles.data.processing_functions{end+1} = processing_function;
handles.data.data_status = min(2,handles.data.data_status);
load_handles(handles);
guidata(hObject,handles)

% --------------------------------------------------------------------
function add_epochs_filter_Callback(hObject, eventdata, handles)
% hObject    handle to add_epochs_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles.data, 'epoched_data')
    return
end
processing_function.str = 'w_filter_epochs';
processing_function.pos = get_first_coincidence(get(handles.processing_select_menu, 'String'),'w_filter_epochs');
processing_function.params = {strjoin(unique({handles.data.epoched_data.epoch.eventtype}))};
input = handles.data.processing_input{processing_function.pos};
if ~isempty(input(:,1))
    answer = dynamic_inputdlg(input(:,1),'Input', 1, processing_function.params);
    if ~isempty(answer)
        processing_function.params = answer;
    else
        return
    end
end
handles.data.processing_functions{end+1} = processing_function;
handles.data.data_status = min(2,handles.data.data_status);
load_handles(handles);
guidata(hObject,handles)

% --- Executes on selection change in processing_functions_list.
function processing_functions_list_Callback(hObject, eventdata, handles)
handles.data.current_processing_function_to_delete = get(hObject,'Value');
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function processing_functions_list_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in remove_processing_function_button.
function remove_processing_function_button_Callback(hObject, eventdata, handles)
if isfield(handles.data,'current_processing_function_to_delete')
    if ~isempty(handles.data.processing_functions(:))
        handles.data.processing_functions(:, handles.data.current_processing_function_to_delete) = [];
        handles.data.data_status = min(2,handles.data.data_status);
    end
    load_handles(handles);
    guidata(hObject,handles)
end



% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function plot_preprocessing_Callback(hObject, eventdata, handles)
    if handles.data.data_status < 1
        handles = run_preprocessing(hObject, handles);
    end
    figure
    plot(handles.data.preprocessed_data);



% --------------------------------------------------------------------
function plot_processing_Callback(hObject, eventdata, handles)
    if handles.data.data_status < 3
        handles = run_processing(hObject, handles);
    end
    figure
    plot(handles.data.processed_data);

% --------------------------------------------------------------------
function plot_scroll_epoching_Callback(hObject, eventdata, handles)
    if handles.data.data_status < 2
        handles = run_epoching(hObject, handles);
    end
    if isfield(handles.data,'epoched_data')
        pop_eegplot( handles.data.epoched_data, 1, 1, 1);
    end


% --- Executes on key press with focus on epoching_variables and none of its controls.
function epoching_variables_KeyPressFcn(hObject, eventdata, handles)
handles.data.data_status = min(1,handles.data.data_status);
load_handles(handles);


% --------------------------------------------------------------------
function run_processing_Callback(hObject, eventdata, handles)
    handles = run_processing(hObject, handles);
    if handles.data.data_status == 3
        load_handles(handles);
        msgbox('Processing finished!','Success')
    end

% --------------------------------------------------------------------
function run_epoching_Callback(hObject, eventdata, handles)
    handles = run_epoching(hObject, handles);
    if handles.data.data_status == 2
        load_handles(handles);
        msgbox('Epoching finished!','Success')
    end

% --------------------------------------------------------------------
function run_preprocessing_Callback(hObject, eventdata, handles)
    handles = run_preprocessing(hObject, handles);
    if handles.data.data_status == 1
        load_handles(handles);
        msgbox('Preprocessing finished!','Success')
    end

function [handles] = run_preprocessing(hObject, handles)
    if ~isfield(handles.data,'path')
        [file_name, path] = uiputfile('*.mat','Save Workspace As');
        if file_name
            handles.data.file_name = file_name;
            handles.data.path = path;
            save_file( handles);
            load_handles(handles);
            guidata(hObject,handles)
        else
            msgbox('You must save the project before executing scripts.','Error')
            return
        end
    end
    [result, data] = execute_preprocessing(handles);
    handles.data = data;
    handles.data.preprocessed_data = result;
    handles.data.data_status = 1;
    guidata(hObject,handles)

function [handles] = run_epoching(hObject, handles)
    if handles.data.data_status < 1
        handles = run_preprocessing(hObject, handles);
    end
    [result, data] = execute_epoching(handles);
    handles.data = data;
    handles.data.epoched_data = result;
    handles.data.data_status = 2;
    guidata(hObject,handles)

function [handles] = run_processing(hObject, handles)
    if handles.data.data_status < 2
        handles = run_epoching(hObject, handles);
    end
    [result, data] = execute_processing(handles);
    handles.data = data;
    handles.data.processed_data = result;
    handles.data.data_status = 3;
    guidata(hObject,handles)


% --------------------------------------------------------------------
function load_project_Callback(hObject, eventdata, handles)
[file_name, path] = uigetfile;
if file_name
    handles = load_file( file_name, path, handles);
    load_handles(handles);
    guidata(hObject,handles)
end


% --------------------------------------------------------------------
function plot_scroll_Callback(hObject, eventdata, handles)
if isfield(handles.data,'EEG')
    pop_eegplot( handles.data.EEG, 1, 1, 1);
else
    msgbox('No .set file has been loaded.','Warning')
end

% --------------------------------------------------------------------
function current_plot_scroll_Callback(hObject, eventdata, handles)
switch handles.data.data_status
    case 0
        if isfield(handles.data,'EEG')
            pop_eegplot( handles.data.EEG, 1, 1, 1);
        else
            msgbox('No .set file has been loaded.','Warning')
        end
    case 1
        pop_eegplot( handles.data.preprocessed_data, 1, 1, 1);
    case 2
        pop_eegplot( handles.data.epoched_data, 1, 1, 1);
    case 3
        pop_eegplot( handles.data.processed_data, 1, 1, 1);
end

% --------------------------------------------------------------------
function load_set_Callback(hObject, eventdata, handles)
[file_name, file_path] = uigetfile('*.set','Select the .set file');
if file_name
    EEG = pop_loadset('filename',file_name,'filepath', file_path);
    EEG = eeg_checkset( EEG );
    handles.data.EEG = EEG;
    handles.data.set_file_name = file_name;
    handles.data.data_status = 0;
    load_handles(handles);
    guidata(hObject,handles)
end


function [data] = calculate_channels_to_discard(handles)
    if isfield(handles.data,'preprocessed_data')
        eeg_data = handles.data.preprocessed_data.data;
    else
        eeg_data = handles.data.EEG.data;
    end
    addpath(fullfile(handles.data.ieeglab_path, 'preprocessing'));
    [channels_to_discard, median_variance, jumps, nr_jumps]  = get_channels_to_discard(eeg_data, 200);
    data.channels_to_discard = channels_to_discard;
    data.median_variance = median_variance;
    data.jumps = jumps;
    data.nr_jumps = nr_jumps;
    disp(channels_to_discard);
    assignin('base', 'channels_to_discard', data)
    

% --------------------------------------------------------------------
function calculate_channels_to_discard_Callback(hObject, eventdata, handles)
calculate_channels_to_discard(handles);

% --------------------------------------------------------------------
function set_channels_to_discard_Callback(hObject, eventdata, handles)
if ~ isfield(handles.data,'channels_to_discard')
    data = calculate_channels_to_discard(handles);
    handles.data.channels_to_discard = data.channels_to_discard;
end
answer = inputdlg('Channels to discard','Input', 5, {mat2str(handles.data.channels_to_discard)});
if ~isempty(answer)
    handles.data.channels_to_discard = str2num(answer{1});
    handles.data.data_status = 0;
    load_handles(handles);
    msgbox('Channels will be discarded after running w_preprocess.','Success')
end
guidata(hObject,handles)


% --------------------------------------------------------------------
function plot_scroll_preprocessed_Callback(hObject, eventdata, handles)
if isfield(handles.data,'preprocessed_data')
    pop_eegplot( handles.data.preprocessed_data, 1, 1, 1);
end


% --------------------------------------------------------------------
function channels_spectra_and_map_Callback(hObject, eventdata, handles)
if isfield(handles.data,'preprocessed_data')
    EEG = handles.data.preprocessed_data;
    answer = inputdlg({'percent data to sample', 'plotting frequency range'},'Input', 1, {'15', '[2 25]'});
    if ~isempty(answer)
        figure; pop_spectopo(EEG, 1, [EEG.xmin*1000 EEG.xmax*1000], 'EEG' , 'percent', str2num(answer{1}), 'freqrange',str2num(answer{2}),'electrodes','off');
    end
else
    msgbox('You should run preprocessing before.','Warning')
end

function close_GUI(hObject, eventdata, handles)
    handles = guidata(hObject);
    selection = questdlg('Do you want to save changes?', ...
	'Warning', ...
	'Yes','No','No');
    switch selection
        case 'Yes'
            if ~isfield(handles.data,'path')
                [file_name, path] = uiputfile;
                if file_name
                    handles.data.file_name = file_name;
                    handles.data.path = path;
                else
                    return
                end
            end
            save_file(handles);
            delete(gcf)
        case 'No'
            delete(gcf)
    end


% --------------------------------------------------------------------
function load_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function view_log_Callback(hObject, eventdata, handles)
    if isfield(handles.data,'path')
        display_text_file(handles.data.path);
    else
        msgbox('Nothing has been logged yet.','Warning')
    end


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function load_edf_Callback(hObject, eventdata, handles)
%     [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%     EEG = pop_biosig('/Users/maximo/Downloads/Florencia_Ignacio_24082014_EFP_2.edf');
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
    [file_name, file_path] = uigetfile('*.edf','Select the .edf file');
    if file_name
        EEG = pop_biosig(fullfile(file_path, file_name));
        EEG = eeg_checkset( EEG );
        handles.data.EEG = EEG;
        handles.data.set_file_name = file_name;
        handles.data.data_status = 0;
        load_handles(handles);
        guidata(hObject,handles)
    end
