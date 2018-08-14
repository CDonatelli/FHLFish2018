function [time, errors] = wobbleWaveKinematicsBatch()
    time = [];
    errors = struct;
    errors.start = 'empty';
    list = dir('*.mat');
    for i = 1:length(list)
           t = [];  tic
       try
           NameStr = list(i).name;
           Struct = load(NameStr);
           Struct = Struct.Struct;
%            if isfield(Struct,'wobblwAmps') == 1
%                 disp([NameStr, ' has already been analyzed.']);
%            else
                disp(NameStr);
                Struct = addBodyAmplitude(Struct);
                save(NameStr, 'Struct');
                t = toc;
%            end
  
       catch err
           errors.(NameStr(1:8)) = getReport(err);
           disp([NameStr, ' had an error']);
       end
       close all
           time = [time; t];
    end

end
