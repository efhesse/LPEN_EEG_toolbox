function [ str, input_variables ] = read_file_by_lines( filename )
%READ_FILE_BY_LINES Reads a multilined file and returns a string with the
%file lines.
%   Reads a multilined file and returns a string with the
%   file lines.
    fileID = fopen(filename,'r');
    str = '';
    i = 1;
    tline = fgetl(fileID);
    input_variables = cell(0);
    while ischar(tline)
        split_tline = strsplit(tline, '| ');
        if isequal(str,'')
            str = split_tline(1);
        else
            str{length(str)+1} = split_tline{1};
        end
        aux = cell(length(split_tline)-1, 2);
        for j = 1:length(split_tline)-1
            var = strsplit(split_tline{j+1}, '= ');
            aux{j,1} = var{1};
            if length(var) == 2
                aux{j,2} = var{2};
            else
                aux{j,2} = '';
            end
        end
        input_variables{i} = aux;
        tline = fgetl(fileID);
        i = i + 1;
    end
    fclose(fileID);
end

