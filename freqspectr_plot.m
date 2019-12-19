directory = biopil_rdir(fullfile('D:\', 'CohortI_c1project', '4_Raw_data', 'ASSR', '**', '*.bdf*')); %bdf or bdf.gz % ** D vervangen afhankelijk van naam schijf
if isempty(directory)
    error('No BDF files found, please check directory name.')
end

params = { ...
    'StimulationSampleRate', 32000, ...
    'StimulationFramesPerEpoch', 32768, ...
    'ReferenceChannels', 'Cz', ... % 'all'
    'Channels', 'all'} ;
% meestal Cz als referentie tegenover alle anderen

%% Read files
for k = 1:numel(directory)
    fprintf('Processing %s...\n', directory(k).name)
    try        
        bdfname = [directory(k).name];
        barename = regexprep(bdfname, '\.bdf$', '');
        [filepath, filename, fileext] = fileparts(barename);
        filenameparts = biopil_strsplit(filename, '_');
        
        if strfind(filename, 'QUIET')
            fprintf('Skipping because QUIET condition\n')
            continue
        end
        
        freq = str2double(regexp(filenameparts(2), '\d+', 'match', 'once'));
        plotname = [barename '_freqspectr' '.png'];
        [filedir, ~, ~] = fileparts(barename);
        [~, lastdir, ~] = fileparts(filedir);
        
%          if exist(csvname, 'file')
%              fprintf('Skipping because csv file exists, remove to rerun\n')
%              continue
%          end
        
        clear s
       
        s = biopil_raw_data( ...
            'FileName', bdfname, ...
            'MultipleEpochs', 'report', ... % 'report' for weird measurements
            'IrregularConditions', 'report', ...
            'Channels', 'all', ...
            'TriggerOffsets', false);
        s.FileHeader.Subject = lastdir;
        
        %% plot frequency spectrum
        s = biopil_plotfd(s, ...
            'ReferenceChannel', 'Cz');

      catch ME
        fprintf('Error during processing of %s:\n', directory(k).name)
        disp(ME)
    end
    
    close all
end
fprintf('Completed\n')      
        
        