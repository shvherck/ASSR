directory = biopil_rdir(fullfile('K:\GBW-0071_DYSCO_Project', 'CohortI_c1project', '4_Raw_data', 'ASSR', '**', '*.bdf*')); %bdf or bdf.gz % ** vervangen door pre of post afhankelijk van wat je wil
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
        csvname =[barename '.csv'];
        plotname = [barename '.png'];
        [filedir, ~, ~] = fileparts(barename);
        [~, lastdir, ~] = fileparts(filedir);
        
        if exist(csvname, 'file')
            fprintf('Skipping because csv file exists, remove to rerun\n')
            continue
        end
        
        clear s
        
        fprintf('Analyzing at %d Hz\n', freq)
        
        s = biopil_raw_data( ...
            'FileName', bdfname, ...
            'MultipleEpochs', 'report', ... % 'report' for weird measurements
            'IrregularConditions', 'report', ...
            'Channels', 'all', ...
            'TriggerOffsets', false);
        s.FileHeader.Subject = lastdir;
        
        % Remove spurious triggers introduced by jitter on the falling edge
        % of the real trigger pulses
        spurioustriggers = find(diff(s.RawData.Triggers) < 100) + 1;
        if numel(spurioustriggers) > 0
            fprintf('Removing %d spurious triggers!\n', numel(spurioustriggers))
            s.RawData.Triggers(spurioustriggers) = [];
            s.RawData.TriggerNames(spurioustriggers) = [];
        end
        
        s = biopil_filter(s, ...
            'HighPass', 2, ...
            'LowPass', []);
        
        s = biopil_triggersclean(s); % remove corrupted triggers
        
        s = biopil_epochs(s);
        s.RawData.EegData = [];
        
        disp(s.Epochs.FramesPerEpoch); % Has to be 8388(triggerlengte)
        disp(s.Epochs.EpochCount);
        s = biopil_epochs_reject(s, ...
            'MaximumEpochCount', 448, ...
            'Percentile', 90); %,'BadChannels', {'Fp1'});
        
        %% Statistics
        % s = biopil_assr_f(s, ...
        %     'EpochsPerSweep', 64, ...
        %     'NoiseBins', 120, ... % number of noise bins for f test
        %     bdfparams);
        
        s = biopil_assr_ht2(s, ...
            params, ...
            'StimulationFrequencies', freq + [-1 0 1] ...
            );
        
        %% create csv file and plot
        s = biopil_csvexport(s);
        
        s = biopil_csvplot(s, ...
            'DetailPlots', {'statistic', 'snr', 'response', 'noise', 'phase'}, ...
            'FigureFileName', plotname);%cvsname
        
        biopil_csvwrite(s, ...
            'FileName', csvname);
        
    catch ME
        fprintf('Error during processing of %s:\n', directory(k).name)
        disp(ME)
    end
    
    close all
end
fprintf('Completed\n')