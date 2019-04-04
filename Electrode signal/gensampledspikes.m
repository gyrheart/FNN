function [retsamplevalues, actualpeaktimes] = gensampledspikes(spiketimes, spiketemplate, varargin) ;
% takes in anarray of spike times (1-d array), in seconds, and an array describing a
% spike template (2d array, N_items * 2), where (:,1) is the times (in ms)
% and (:,2) is the membrane voltage, and produces from this a 1-d array of
% the values of the membrane voltage at sample times. The sampling rate
% defaults to 100,000 (10 microseconds), but can be set.
%
% LSS started 9 June 2005.
% LSS 24 July 2017: in call to interp1, cubic changed to pchip
%

% Read in arguments
SampleRate = 100000 ;
RestValue = spiketemplate(1,2) ; % default is first value in spike template
MinVal = -0.1 ;
MaxVal = 1.0 ;
% default end time final maximum length is last spike time + total duration of the
% spiketemplate, rounded up to a spmple point

if ~isempty(spiketimes)
    endtime = spiketimes(length(spiketimes)) + spiketemplate(length(spiketemplate), 1)/1000 ;
    returnendtime = endtime ;
end
for i=1:(nargin-2)
    if ischar(varargin{i})
        switch varargin{i}
            case 'SampleRate', SampleRate = varargin{i+1};
            case 'RestValue', RestValue = varargin{i+1};
            case 'MinVal', MinVal = varargin{i+1};
            case 'MaxVal', MaxVal = varargin{i+1};   
            case 'EndTime' 
                endtime = varargin{i+1} + spiketemplate(length(spiketemplate), 1)/1000 ; 
                returnendtime = varargin{i+1} ;
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
sampletemplate = interp1(spiketemplate(:,1)/1000, spiketemplate(:,2), ...
    spiketemplate(1,1)/1000:1/SampleRate:spiketemplate(length(spiketemplate), 1)/1000, 'pchip') ;
templatelength = length(sampletemplate) ;

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
