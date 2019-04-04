function [finalsignal actualpeaks stimes_1 rngreturnstatus] = makenoisysamples(varargin)
%
% this functionputs together all the functions in this directory, to
% produce a noisy spike train. It really *ought* to be a GUI, but well...it
% isn't (yet).
%
% Instead, all the parameters that *would* be filled in in the GUI are here at the start of the function.
%
% returns (i) the noisy signal aand (ii) the actual (official!) spike
% times. These are the peaks of the spikes. (iii) The spike times used to generate the spike train. 
% Note that the spike times are the *starting times* of the spike
% templates, so that the actual peak will be some little time later. This
% last set are included so that one can re-do data sets with the same
% primary spikes, but different noise levels etc. See varargin below.
%
% varargin: just 2 so far: ReuseTimes, followed by the actual times of
% spike starts (so that one can generate a number of sets of spikes with
% different noise levels etc., but the *same* underlying sets of primary
% spikes
%  ReuseRNGstate: re use the same random number generator state
%
% LSS started 8 June 2005
% First release 14 June 2005.
% current release 21 Oct 2005.
%
sample_rate = 100000 ;

templatename = 'spike_1n.csv' ; % template for spike for primary neuron
duration = 0.05 ; % length of the signal: NB: approximate. actual length may be a little longer
spikedist = 'P' ; % 'P' for poission, 'G' for Gaussian
% for poission distribution only
poissionnumber = 6 ;
poissmeanISI = 0.01 ;
%for gaussian only
gaussmeanITD = 0.01 ;
gaussstdev = 0.002 ;
% for differential
dslength = 15 ; % used in smoothing during diferentiation. Is a number of samples.

% mixing coefficients for the primary neural signal.
origmixcoeffs_orig = 0 ; % plain (intracellular) signal strength
origmixcoeffs_d = -0.1 ; % differentiated signal strength
origmixcoeffs_dd = 1 ; % twice differentiated signal strength

% for jittered noise
n_jitter = 3; ; %number of jittered versions of the original spike train
jitteroverallsize_orig = 0.1 ; %overall size of the original jittered signal
jitteroverallsize_d = 0.1 ; %overall size of the differentiated jittered signal
jitteroverallsize_dd = 0.1 ; %overall size of the twice differentiated jittered signal
% The arrays below must ne at least n_jittel long
jprobspikes = [0.99 0.98 0.97] ; % probability that an inpout spike results in an output spike
jjstd = [0.004 0.003 0.005] ; % standard deviation of the jitter (mean = 0)
% jittered noise mixing coefficients
jmixcoeffs_orig = [0 0 0] * jitteroverallsize_orig  ; % plain (intracellular signal strength
jmixcoeffs_d = [-0.1 -0.095 -0.112] * jitteroverallsize_d ; % for 1st differential
jmixcoeffs_dd = [1 1.1 0.9] * jitteroverallsize_dd; % for 2nd differential

% for uncorrelated spike noise
n_uncorr = 2 ; % number of uncorrelated spike trains
u_overallsize_orig = 0.1 ; %overall size of the original un-jittered signal
u_overallsize_d = 0.1 ; %overall size of the differentiated un-jittered signal
u_overallsize_dd = 0.1 ; %overall size of the twice differentiated un-jittered signal
% All arrays below must be at least n_uncorr long
u_templatename = ['spike_1n.csv'; 'spike_1n.csv'] ; % tempate for these spikes (N.B. not really implemented yet)
u_spikedist = ['P' 'P'] ; % distribution type for these spike trains
u_poissionnumber = [4 7] ; % number for use in generating Poisson distribution
u_poissmeanISI = [0.01 0.015] ; % Mean ISI for each poisson distribution
u_gaussmeanITD = [0.01 0.01] ; % Mean Interspike time difference (shoud be ISI!)
u_gaussstdev = [0.005 0.005] ; % STD thereof
% uncorrelated spikes mixing coefficients
umixcoeffs_orig = [0 0] * u_overallsize_orig ; % for intracellular signal
umixcoeffs_d = [-0.1 -0.15] * u_overallsize_d ; % for 1st differential
umixcoeffs_dd = [0.9 1.1] * u_overallsize_dd ; % for 2nd differential

% noise level
% in dB (for awgn function)
noise_snr = 20 ;

% final output range
finalminval = -1.0 ;
finalmaxval = 1.0 ;

% end of parameters!


% read in the spike template
spiketemplate = csvread(templatename) ;

% generate some spike times
if spikedist == 'P'
    stimes_1 = genspikespoisson(duration, poissionnumber, poissmeanISI) ;
else if spikedist == 'G'
        stimes_1 = genspikesgaussion(duration, gaussmeanITD, gaussstdev) ;
    end
end

% overwrite stimes_1 if ReuseTimes is provided
for i=1:(nargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'ReuseTimes', stimes_1 = varargin{i+1};
            case 'ReuseRNGstate', rngstate = varargin{i+1};
        end
    end
end

% reset the random number generator if appropriate
if exist('rngstate') == 1
    rand('state',rngstate) ;
else
    rand('state',sum(100*clock)) ;
end
rngreturnstatus = rand('state') ; % in case we might want to use it again


% generate some sampled spike data
[ssamples_1 actualpeaks] = gensampledspikes(stimes_1, spiketemplate, 'EndTime', duration, 'SampleRate', sample_rate) ;

% generate the differential and the second differential of both of these.
ssamples_1d = simplediff(ssamples_1, dslength) ;
ssamples_1dd = simplediff(ssamples_1d, dslength) ;
% and normalise them
ssamples_1n = setlimits(ssamples_1, -0.5, 0.5) ;
ssamples_1dn = setlimits(ssamples_1d, -0.5, 0.5) ;
ssamples_1ddn = setlimits(ssamples_1dd, -0.5, 0.5) ;

% generate some jittered data
for nj = 1: n_jitter
    % generate jittered noisy spike trains
    jstimes = probjitter(stimes_1, jprobspikes(nj), jjstd(nj), duration );
    [jsamples temp] = gensampledspikes(jstimes, spiketemplate, 'EndTime', duration, 'SampleRate', sample_rate) ;
    % and its derivatives
    jsamples_d = simplediff(jsamples, dslength) ;
    jsamples_dd = simplediff(jsamples_d, dslength) ;
    % and normalise the derivatives
    jsamples_n(nj, :) = setlimits(jsamples, -0.5, 0.5) ;
    jsamples_dn(nj, :) = setlimits(jsamples_d, -0.5, 0.5) ;
    jsamples_ddn(nj, :) = setlimits(jsamples_dd, -0.5, 0.5) ;
end

% generate some uncorrelated spikes too
for nu=1:n_uncorr
    if u_spikedist(nu) == 'P'
        u_stimes = genspikespoisson(duration, u_poissionnumber(nu), u_poissmeanISI(nu)) ;
    else if u_spikedist(nu) == 'G'
            u_stimes = genspikesgaussion(duration, u_gaussmeanITD(nu), u_gaussstdev(nu)) ;
        end
    end
    % get in the spike template (might be different: if it is, may need to fix duration up later)
    u_spiketemplate = csvread(u_templatename(nu,:)) ;
    % generate some sampled spike data
    [u_ssamples temp] = gensampledspikes(u_stimes, u_spiketemplate, 'EndTime', duration, 'SampleRate', sample_rate) ;
    u_ssamples_d = simplediff(u_ssamples, dslength) ;
    u_ssamples_dd = simplediff(u_ssamples_d, dslength) ;
    % and normalise them
    u_samples_n(nu,:) = setlimits(u_ssamples, -0.5, 0.5) ;
    u_samples_dn(nu,:) = setlimits(u_ssamples_d, -0.5, 0.5) ;
    u_samples_ddn(nu,:) = setlimits(u_ssamples_dd, -0.5, 0.5) ;
end

% save a little space
clear u_ssamples u_ssamples_d u_ssamples_dd 
% now mix these together
corrspikenoise = zeros([1 length(ssamples_1)] );
for i_corr = 1:n_jitter
    corrspikenoise = corrspikenoise + jmixcoeffs_orig(i_corr) * jsamples_n(i_corr, :) ;
    corrspikenoise = corrspikenoise + jmixcoeffs_d(i_corr) * jsamples_dn(i_corr, :) ;
    corrspikenoise = corrspikenoise + jmixcoeffs_dd(i_corr) * jsamples_ddn(i_corr, :) ;
end
% save some more space
clear jsamples_n jsamples_dn jsamples_ddn

uncorrspikenoise = zeros([1 length(ssamples_1)] );
for i_uncorr=1:n_uncorr
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_orig(i_uncorr) * u_samples_n(i_uncorr,:) ;
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_d(i_uncorr) * u_samples_dn(i_uncorr,:) ;
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_dd(i_uncorr) * u_samples_ddn(i_uncorr,:) ;
end

% clear some space
clear u_samples_n u_samples_dn u_samples_ddn

origsignal = origmixcoeffs_orig * ssamples_1n +  ...
    origmixcoeffs_d * ssamples_1dn + origmixcoeffs_dd * ssamples_1ddn ;

noisysignal = origsignal + corrspikenoise + uncorrspikenoise ; 
% noisysignal = origsignal ;
% add some white noise
almostfinalsignal = awgn(noisysignal, noise_snr, 'measured') ;

% finally, scale the final signal 
finalsignal = setlimits(almostfinalsignal, finalminval, finalmaxval) ;

end

