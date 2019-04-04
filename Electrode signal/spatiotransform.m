function troutput = spatiotransform(samples, samplerate, varargin)
% essentially convolves the DelayCoefficients withe ssample input to
% produce the output. However, the samples are at the sample rate, but the
% Delay Coefficients are at a pitch of DelayDelta.
%
% LSS started 22 2 2006
% LSS 24 July 2017: in call to interp1, cubic changed to pchip
%
% default values: used really to make tersting simpler. Normally
% overwritten
DelayDelta = 3.0/samplerate ; 
DelayNumber = 10 ;
DelayCoefficients = [0.1 0.2 0.3 0.4 0.5 0.4 0.3 0.2 0.1 0.05] ;
for i=1:(nargin-2)
    if ischar(varargin{i})
        switch varargin{i}
            case 'DelayDelta', DelayDelta =  varargin{i+1};
            case 'DelayNumber', DelayNumber = varargin{i+1};
            case 'DelayCoefficients', DelayCoefficients = varargin{i+1};
        end
    end
end
%
% Form the DelayCoefficients into something with the same pitch (i.e.
% sampling rate) as the input samples
dcoeffs_length = ceil((DelayDelta * DelayNumber  * samplerate) + 1);
dcoeffs = interp1([0:DelayNumber-1]*DelayDelta, DelayCoefficients, [0:dcoeffs_length-1]/samplerate, 'pchip', 0) ;
%
troutput = conv(dcoeffs, samples) ;

%
