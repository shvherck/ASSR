function s = eegorl_AvgData(s, varargin)
%EEGORL_AVGDATA Combine channels in time domain
%
%  Combine channels in time domain. Overwrite channel with
%  averaged channels
%
%   Input:  biopil struct
%
%   Output: biopil struct
%
%   Usage:  s=eegorl_AvgData(s, 'CombineSet', {'P8','P9'})
%
%   Author: Robert Luke


%% Check input data

p = inputParser;
p.KeepUnmatched = true;
p.StructExpand  = true;
p.addParamValue('CombineSet'       , {'Sophie'}              , @iscell   );
p.addParamValue('OverwriteChannel' , 'Fpz'                   , @isstr    );
p.addParamValue('NewName'          , 'AvgData'               , @isstr    );
p.addParamValue('BadChannels'      , {'Status'}              , @iscell   );
p.addParamValue('Verbose'          , 1                       , @isnumeric);
p.parse( varargin{:})
params = p.Results;


%% Check electrodes to average and overwrite

if strcmpi(params.CombineSet, 'Sophie')
    selected_electrodes_all = {'TP7', 'P9',  'P7', 'P5', 'P3', 'P1', 'PO7', 'PO3', 'O1', ...
                               'TP8', 'P10', 'P8', 'P6', 'P4', 'P2', 'PO8', 'PO4', 'O2'};
elseif strcmpi(params.CombineSet, 'CIset')
    selected_electrodes_all = {'P9', 'P10', 'PO7', 'PO8', 'PO3', 'PO4', 'CP5', 'CP6', 'O1', 'O2', 'Oz', 'Iz', 'Pz'};
elseif strcmpi(params.CombineSet, 'SophieLeft')
    selected_electrodes_all = {'TP7', 'P9',  'P7', 'P5', 'P3', 'P1', 'PO7', 'PO3', 'O1'};
elseif strcmpi(params.CombineSet, 'SophieRight')
    selected_electrodes_all = {'TP8', 'P10', 'P8', 'P6', 'P4', 'P2', 'PO8', 'PO4', 'O2'};
elseif strcmpi(params.CombineSet, 'Frontal')
    selected_electrodes_all = {'Fp1', 'Fpz', 'Fp2', 'AF1', 'AF3', 'AFz', 'AF4', 'AF8'};
else
    % Check can find all the channels
    num_valid_channels  = sum(ismember({s.FileHeader.Channels.Label}, params.CombineSet));
    if num_valid_channels == numel(params.CombineSet)
        % All channels are valid
        selected_electrodes_all = params.CombineSet;
    else
        warning('Some of the channels you wish to combine cant be found. Aborting');
        return
    end
end

% Remove bad electrodes from set
selected_electrodes_all(ismember(selected_electrodes_all, params.BadChannels)) = [];

% Find electrode indicies
selected_electrodes  = ismember(lower({s.FileHeader.Channels.Label}), lower(selected_electrodes_all));
overwrite_channel    = ismember(lower({s.FileHeader.Channels.Label}), lower(params.OverwriteChannel));

% Check only one overwrite channel
if sum(overwrite_channel)~=1
    warning('Error with channel to over write. Aborting');
    return
end

% User feedback if requested
if params.Verbose
    fprintf('Combining %2d channels and replacing channel %s=%d\n', ...
        numel(selected_electrodes_all), params.OverwriteChannel, find(overwrite_channel));
end
if isfield(s, 'Epochs')
    warning('You already have epochs calculated, these will no include the averaged data');
end


%% Combine channels

s.RawData.EegData(:,overwrite_channel) = mean(s.RawData.EegData(:,selected_electrodes),2);
s.FileHeader.Channels(overwrite_channel).Label = params.NewName;


%% Save data in biopil structure

% This will only work once

s.AvgData = params;
s.AvgData.OverwriteChannelIndex = overwrite_channel;
s.AvgData.CombineSetIndex       = selected_electrodes;
