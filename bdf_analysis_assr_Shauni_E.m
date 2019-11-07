%% Create directory and define parameters
directory = biopil_rdir(fullfile('D:\', 'CohortI_c1project', '4_Raw_data', 'ASSR', '**', '*.bdf*')); %bdf or bdf.gz %hard drive may be E: or D: (check on pc)
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
        %% Define parameters
        bdfname = [directory(k).name];
        barename = regexprep(bdfname, '\.bdf$', '');
        [filepath, filename, fileext] = fileparts(barename); %fileparts returns filepath, name & extension
        filenameparts = biopil_strsplit(filename, '_');
        
        if strfind(filename, 'QUIET')
            fprintf('Skipping because QUIET condition\n')
            continue
        end
        
        freq = str2double(regexp(filenameparts(2), '\d+', 'match', 'once'));
        csvname =[barename '.csv'];
        plotname = [barename '.png'];
        freqname = [barename '_freqspectr' '.png'];
        [filedir, ~, ~] = fileparts(barename);
        [~, lastdir, ~] = fileparts(filedir);
        
        if exist(csvname, 'file')
            fprintf('Skipping because csv file exists, remove to rerun\n')
            continue
        end
        
        clear s
        
        fprintf('Analyzing at %d Hz\n', freq)
       
        %% Read raw data
        s = biopil_raw_data( ...
            'FileName', bdfname, ...
            'MultipleEpochs', 'report', ... % 'report' for weird measurements
            'IrregularConditions', 'report', ...
            'Channels', 'all', ...
            'TriggerOffsets', false);
        s.FileHeader.Subject = lastdir;
        
        %% Clean data
        % Remove spurious triggers introduced by jitter on the falling edge
        % of the real trigger pulses
        spurioustriggers = find(diff(s.RawData.Triggers) < 100) + 1;
        if numel(spurioustriggers) > 0
            fprintf('Removing %d spurious triggers!\n', numel(spurioustriggers))
            s.RawData.Triggers(spurioustriggers) = [];
            s.RawData.TriggerNames(spurioustriggers) = [];
        end
        
        s=eegorl_AvgData(s, 'CombineSet', {'SophieLeft'},'OverwriteChannel' , 'Fpz','NewName', 'Left');
        s=eegorl_AvgData(s, 'CombineSet', {'SophieRight'},'OverwriteChannel' , 'AFz','NewName', 'Right');
        s=eegorl_AvgData(s, 'CombineSet', {'Sophie'},'OverwriteChannel' , 'Fz', 'NewName', 'All');

        s = biopil_filter(s, ...
            'HighPass', 2, ...
            'LowPass', []);
        
        s = biopil_triggersclean(s); % remove corrupted triggers
        
        s = biopil_epochs(s);
        s.RawData.EegData = [];
       
        %remove bad recording  (see baddata_removal_efficient.m for remarks
        %on working script
        data = readtable('baddata.xlsx');
        str = string(data{:,1});                
        pattern = string(filename);             
        matches = strfind(str,pattern);         
        tf = any(vertcat(matches{:}));                                         
        if tf == 1                                   
            index = string(data{:,1}) == filename;   
            timepoint1 = double(data{index,2}); %starting point of bad recording     
            timepoint2 = double(data{index,3}); %finishing point of bad recording    
            fprintf(filename);
            disp(timepoint1);
            disp(timepoint2);
            sample_1 = timepoint1*s.FileHeader.SampleRate; %multiply seconds by sampling rate to get the sample at a timepoint
            sample_2 = timepoint2*s.FileHeader.SampleRate;        
            h = floor(sample_1/s.Epochs.FramesPerEpoch); %divide sample by frames per epoch to get the epoch at a timepoint
            hi = ceil(sample_2/s.Epochs.FramesPerEpoch);
            s.Epochs.EegData(:,h:hi,:)=[]; %remove everything between the two epochs
            s.Epochs.EpochCount = length(s.Epochs.EegData(1,:,1)); %change epoch count to the number of remaining epochs
        end

        disp(s.Epochs.FramesPerEpoch); % Has to be 16776(triggerlengte)
        disp(s.Epochs.EpochCount);
        s = biopil_epochs_reject(s, ...
            'MaximumEpochCount', 448, ... % selects the 1st n epochs
            'Percentile', 90); %,'BadChannels', {'Fp1'}); % rejects epochs based on peak-to-peak amplitude
        
        s.RawData.EegData = [];
        s.AvgData = [];

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
        
        %% plot frequency spectrum
        s = biopil_plotfd(s, ...
            'StimulationSampleRate', 32000, ...
            'StimulationFramesPerEpoch', 32768, ...
            'StimulationFrequencies', freq, ...
            'Channels', {'Left','Right','All'}, ...
            'FigureFileName', freqname);

    catch ME
        fprintf('Error during processing of %s:\n', directory(k).name)
        disp(ME)
    end
    
    close all
end
fprintf('Completed\n')