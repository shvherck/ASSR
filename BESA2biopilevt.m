function s = BESA2biopilevt(s,evt)
% Function used to apply epoch rejection performed in BESA to data
% structure in biopil format

% evt - 3D matrix read in to matlab by BESA reader, 1st column contains time stamps of accepted triggers in us
% s   - bioipl struct with EEG data
% 2 Dec 2019 Anna Samsel

 time       = 0.0:1.0/s.FileHeader.SampleRate:((size(s.RawData.EegData,1)-1)/s.FileHeader.SampleRate); % creating time array
 time       = time * 10^(6); % in order to have the same scale as BESA data
 biopil_evt = time(s.RawData.Triggers);% taking timestamps of consequtive biopil triggers

 [val,pos]  = intersect(int32(biopil_evt), evt(:,1));% finding value and position of triggers in biopil and BESA taht occur at the same time
 
% sanity check up
if length(pos) ~= size(evt,1) %because in the file that I got from you this was not the case, but i can imagine due to rounding this might happen in some cases
    error('Length of detected common triggers for biopil and BESA and evt BESA file is different. Please check if time stamps for BESA and biopil match.')
end

% ctreating mask to exclude some triggers from biopil
mask      = logical(zeros(length(s.RawData.Triggers),1));
mask(pos) = logical(1);

% apllying mask to biopil triggers definition, which will get rid of the
% triggers that were removed in BESA analysis pipeline

s.RawData.Triggers = s.RawData.Triggers(mask);
s.RawData.TriggerNames = s.RawData.TriggerNames(mask);
s.FileHeader.RecordCount = length(pos);


end