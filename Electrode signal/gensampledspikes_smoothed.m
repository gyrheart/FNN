function [retsamplevalues actualpeaktimes] = gensampledspikes_smoothed(spiketimes, spiketemplate, varargin) ;
% takes in anarray of spike times (1-d array), in seconds, and an array describing a
% spike template (2d array, N_items * 2), where (:,1) is the times (in ms)
% and (:,2) is the membrane voltage, and produces from this a 1-d array of
% the values of the membrane voltage at sample times. The sampling rate
% defaults to 100,000 (10 microseconds), but can be set.
%
% LSS started 9 June 2005.
%
% new version: Allows multiple weighted echoes of the intracellular potential, emulating
% the effect of the spatial extent of the neuron. If the weights are all the same,
% this models signals from distant neurons for which the actual shape 
% of the neuron is not important, but the duation over which the spike is emitted 
% by different parts of the neuron simply smooths the spikes because the transfer
% function for the whole spatial extent is much the same. On the other
% hand, we can make the weights vary, emulating (e.g.) the neuron being
% near to the electrode, with some parts much nearer than others, so that
% the transfer function differs from part to part of the neuron. Note that
% the actual overall magnitude is ignored, since the signal is scaled at
% the end. 
%
% new version started 26 Jan 2006.

% Read in arguments
SampleRate = 100000 ;
RestValue = 0 ; % default is 0. First value in smoothed spike template is 0.
MinVal = -0.1 ;
MaxVal = 1.0 ;
% default end time final maximum length is last spike time + total duration of the
% spiketemplate, rounded up to a spmple point
% arguments for delay based sommothing
DelayDelta = 0.1 ; % in milliseconds
DelayNumber = 20 ; % number of delayed additions
DelayCoefficients = ones([1 DelayNumber]) ;
%
if ~isempty(spiketimes)
    % include a little extra on the end
    endtime = spiketimes(length(spiketimes)) + spiketemplate(length(spiketemplate), 1)/1000 + 0.02 ;
    returnendtime = endtime ;
end
for i=1:(nargin-2)
    if ischar(varargin{i})
        switch varargin{i}
            case 'SampleRate', SampleRate = varargin{i+1};
            case 'RestValue', RestValue = varargin{i+1};
            case 'MinVal', MinVal = varargin{i+1};
            case 'MaxVal', MaxVal = varargin{i+1};   
            case 'EndTime', 
                endtime = varargin{i+1} + spiketemplate(length(spiketemplate), 1)/1000 ; 
                returnendtime = varargin{i+1} ;
            case 'DelayDelta', DelayDelta =  varargin{i+1};
            case 'DelayNumber', DelayNumber = varargin{i+1};
            case 'DelayCoefficients', DelayCoefficients = varargin{i+1};
        end
    end
end
    
% set up samplevalues array

lastsample = ceil(endtime * SampleRate) ;
returnlastsample = ceil(returnendtime * SampleRate) ; % for actual returned vector

samplevalues = ones([1 lastsample]) * RestValue ; % resting value

if isempty(spiketimes) % simply return if there's no spikes
    actualpeaktimes = [] ;
    return ;
end
% make the spiketemplate into a vector with a values for each sample
templateduration = spiketemplate(length(spiketemplate), 1) - spiketemplate (1,1) ;


% spikepoint determines where (relative to the spike time) spikes get added
% in
if spiketemplate (1,1) == 0
    spikepoint = 0 ;
else
    spikepoint = floor((-spiketemplate(1,1)/templateduration)) * SampleRate ; % but a better idea is to normalise so start is 0.
    disp(['Warning: (gensampledspikes.m) spike generation template does not start at 0']) ;
end
% Use interp1 to generate all the interpolated values
s_template = interp1(spiketemplate(:,1)/1000, spiketemplate(:,2), ...
    spiketemplate(1,1)/1000:1/SampleRate:spiketemplate(length(spiketemplate), 1)/1000, 'cubic') ;
s_templatelength = length(s_template) ;

% normalise this template to be 0 start
s_template = s_template - s_template (1) ; 

% calculate effect of delay
DelaySamples = floor((DelayDelta * SampleRate) /1000) ;
TotalDelaySamples = DelaySamples * DelayNumber ;
% find expected length of sampletemplate
templatelength = s_templatelength + TotalDelaySamples ;
sampletemplate = zeros([1 templatelength]) ;
% initialise sampletemplate
sampletemplate(1:s_templatelength) = s_template ;
% then add in  (weighted) delay effects
for delaynumber = 1:DelayNumber
    sampletemplate(DelaySamples * delaynumber + 1: DelaySamples * delaynumber + s_templatelength) = ...
      sampletemplate(DelaySamples * delaynumber + 1: DelaySamples * delaynumber + s_templatelength) +  ...
      (s_template  * DelayCoefficients(delaynumber)) ;
end
sampletemplate = sampletemplate + s_template (1) ; % add back in (doesn't really matter, as result is normalised at the end)
% find peak of sampletemplate
[temp, spmaxlocn] = max(sampletemplate) ;
actualpeaktimes = zeros(size(spiketimes)) ; % for recording actual peak locations

% calculate gradient of spike at the very start, for possible use later
spikegrad = (sampletemplate(2) - sampletemplate(1)) ; % per sample interval

% add in a sampled spike for each spike time
for spikeno = 1:length(spiketimes)
    % add this spike in
    spikelocation = ceil(spiketimes(spikeno) * SampleRate) ;
    % do we need to smooth the transition from current potential to that at
    % the start of the spike template?
    if samplevalues(spikelocation - spikepoint - 1) == sampletemplate(1)
        samplevalues(spikelocation - spikepoint: spikelocation - spikepoint + templatelength-1) = sampletemplate ;
        actualpeaktimes(spikeno) = (spikelocation - spikepoint + spmaxlocn)/SampleRate ;
    else % smoothing required
        % there's 2 possibilities: (i) previous spike is still in
        % (downswing!) of spike (still in upswing is probably not allowed
        % due to refractory period)
        % or (ii) previous spike is in AHP
        if samplevalues(spikelocation - spikepoint - 1) > RestValue % still in spike itself
            % still has a very sharp corner!
            startindex = find(sampletemplate >= samplevalues(spikelocation - spikepoint - 1), 1) ;
            if isempty(startindex)
                disp(['gensampledspikes.m: startindex empty: problem: please contact lss']) ;
            end
            samplevalues(spikelocation - spikepoint: spikelocation - spikepoint + (templatelength - startindex)) ...
                = sampletemplate(startindex:templatelength) ;
            actualpeaktimes(spikeno) = (spikelocation - spikepoint + (spmaxlocn - startindex))/SampleRate ;
        else % in AHP
            % adjust up to (just below) the start of the sampletemplate
            % using the gradient
            addedlocations = 0 ;
            location = spikelocation - spikepoint -1 ;
            while samplevalues(location) < sampletemplate(1)
                location = location + 1 ;
                samplevalues(location) = samplevalues(location - 1) + spikegrad ;
                addedlocations = addedlocations + 1 ;
            end
            if (location + templatelength-1 <= length(samplevalues))
                samplevalues(location: location + templatelength-1) = sampletemplate ;
            else % in danger of extending the samplevalues array
                samplevalues(location: length(samplevalues)) = sampletemplate(1:length(samplevalues) - location + 1) ;
            end
            actualpeaktimes(spikeno) = (spikelocation - spikepoint + spmaxlocn + addedlocations)/SampleRate ;
        end
    end

end

% normalise the amplitude of the spikes
% stretch linearly so that all values lie between MinVal and MaxVal
samplevalues = setlimits(samplevalues, MinVal, MaxVal);
retsamplevalues = samplevalues(1:returnlastsample) ;
end
