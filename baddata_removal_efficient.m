% name              baddata_removal_efficient.m  
% creator           Shauni Van Herck & Raùl Granados Barbero
% description       used to remove bad data in if-loop by reading in
%                   summary of bad data, checking if the filename exists
%                   in the summary, and if it does, removing the raw data
%                   in between the 2 timepoints defined in the summary
% version.string    MATLAB R2016b            
% platform          x86_64-w64-mingw32                      
% date created      05/11/2019 

filename = 'i035_PULS4_pre';                                               %filename is already definesd in assr script

data = readtable('baddata.xlsx');
str = string(data{:,1});                                                   %str needs to be string
pattern = string(filename);                                                %pattern needs to be string
matches = strfind(str,pattern);                                            %does filename occur in the 1st column of the data table?
tf = any(vertcat(matches{:}));                                             %returns 1 if filename occurs in data table and 0 if not                               
if tf == 1                                                                 %if filename occurs in the data table
    index = string(data{:,1}) == filename;                                 %find the row for filename
    timepoint1 = double(data{index,2});                                    %define timepoints as values in column 2 & 3 of the row for filename
    timepoint2 = double(data{index,3});                                    %needs to be double for future operations
    fprintf(filename);
    disp(timepoint1);
    disp(timepoint2);
    sample_1 = timepoint1*s.FileHeader.SampleRate;                         %multiply seconds by sampling rate to get the sample at a timepoint
    sample_2 = timepoint2*s.FileHeader.SampleRate;        
    h = floor(sample_1/s.Epochs.FramesPerEpoch);                           %divide sample by frames per epoch to get the epoch at a timepoint
    hi = ceil(sample_2/s.Epochs.FramesPerEpoch);
    s.Epochs.EegData(:,h:hi,:)=[];                                         %remove everything between the two epochs
    s.Epochs.EpochCount = length(s.Epochs.EegData(1,:,1));                 %change epoch count to the number of remaining epochs
end

%% Remarks

% to remove the 1st x seconds of a recording you need timepoint1 2s, because this corresponds to epoch 1 - otherwise script will fail
