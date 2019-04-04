function [finalsignal targets rngreturnstatus] = makenoisysamples_2(varargin)
%
% this functionputs together all the functions in this directory, to
% produce a noisy spike train. It really *ought* to be a GUI, but well...it
% isn't (yet).
%
% Instead, all the parameters that *would* be filled in in the GUI are here
% at the start of the function.
%
% returns (i) the noisy signal aand (ii) the actual (official!) spike
% times. These are the peaks of the spikes. (iii) The spike times used to 
% generate the spike train. 
% Note that the spike times are the *starting times* of the spike
% templates, so that the actual peak will be some little time later. This
% last set are included so that one can re-do data sets with the same
% primary spikes, but different noise levels etc. See varargin below.
%
% varargin: just 2 so far: 
% ReuseTimes, followed by the actual times of
% spike starts (so that one can generate a number of sets of spikes with
% different noise levels etc., but the *same* underlying sets of primary
% spikes
% ReuseRNGstate: re use the same random number generator state
%
% New version started 23Oct 2005: have a number of target neurons, each with
% their own parameters, so that one can do spike sorting as well as spike
% detection. Further, allow the correlated spike trains to be correlated to
% one of the target neurons, and allow them to use whatever original spike
% train is required.
%
% LSS started 8 June 2005
% First release 14 June 2005.
% current release 21 Oct 2005.
%
% defaults:  may be overwritten in varargin section
n_targets =1 ; % number of target neurons
duration = 0.05 ; % length of the signal: NB: approximate. actual length may be a little longer
sample_rate = 100000 ;
dslength = 15 ; % used in smoothing during diferentiation. Is a number of samples.

% overwrite stimes_1 if ReuseTimes is provided
for i=1:(nargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'N_Targets'
                n_targets = varargin{i+1};
            case 'ReuseTargets'
                targets = varargin{i+1};
            case 'Duration'
                duration = varargin{i+1};
            case 'SampleRate'
                sample_rate = varargin{i+1};
            case 'DSLength'
                dslength = varargin{i+1};
            case 'ReuseRNGstate'
                rngstate = varargin{i+1};
        end
    end
end

% templates for spike for target (primary) neuron. Note that all file
% names must be same length
templatenames = ['spike_1n.csv',  
    'spike_2n.csv']; 
[n_templates temp] = size(templatenames) ; 

% target neuron information
% needs to be replicated n_targets times
% below is for n_targets = 3

spikedist = ['P' ,% distribution to be used in generation:'P' for poission, 'G' for Gaussian
    'G',
    'P'] ; 
% for poission distribution only
poissionnumber = [6 0 8] ; % poisson values: only used where there's a 'P' in spikedist
poissmeanISI = [0.01 0 0.01];
t_templateid = [1 2 1] ; % template for this spike train
%for gaussian only
gaussmeanITD = [0 0.01 0]; % gaussian values: only used where there's a 'G' in spike dist
gaussstdev = [0 0.002 0 ] ;
% for differential amounts of each of the signals


% mixing coefficients for the primary neural signal.
origmixcoeffs_orig = [0 0 0]; % plain (intracellular) signal strength
origmixcoeffs_d = [-0.1 0.2 0.1]; % differentiated signal strength
origmixcoeffs_dd = [1 0.9 0.87]; % twice differentiated signal strength

% for jittered noise
n_jitter = 3 ; %number of jittered versions of the original spike train
jitter_orig = [1 2 3] ; % which original signal should be jittered
jitteroverallsize_orig = 0.1 ; %overall size of the original jittered signal
jitteroverallsize_d = 0.1 ; %overall size of the differentiated jittered signal
jitteroverallsize_dd = 0.1 ; %overall size of the twice differentiated jittered signal
% The arrays below must ne at least n_jitter long
j_templateid = [2 1] ; % templates for these spikes 
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
u_templateid = [1 2] ; % templates for these spikes (N.B. not really implemented yet)
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


% read in the spike templates (used for targets and uncorrelated spikes)
for i = 1:n_templates
    spiketemplates(i).template = csvread(templatenames(i,:)) ;
    % adjust so that the times are norm alised to start at 0
    [tlength temp] = size(spiketemplates(i).template) ;
    tadj = ones([tlength 1]) * spiketemplates(i).template(1,1) ;
    spiketemplates(i).template(:,1) = spiketemplates(i).template(:,1) - tadj ;
    spiketemplates(i).templatename = templatenames(i,:) ;
end

% reset the random number generator if appropriate
if exist('rngstate') == 1
    rand('state',rngstate) ;
else
    rand('state',sum(100*clock)) ;
end

rngreturnstatus = rand('state') ; % in case we might want to use it again

% if targets has been supplied as varargin, then use them. Otherwise
% generate the target spike times. 
if ~exist('targets')
    % generate the target spike times
    for i = 1:n_targets
        if spikedist(i) == 'P'
            targettimes = genspikespoisson(duration, poissionnumber(i), poissmeanISI(i)) ;
        else if spikedist(i) == 'G'
                targettimes = genspikesgaussian(duration, gaussmeanITD(i), gaussstdev(i)) ;
            end
        end
        targets(i).targettimes = targettimes ;

        % generate some sampled spike data from template t_templateid(i)
        [targets(i).sampletargets targets(i).actualpeaks] = ...
            gensampledspikes(targettimes, spiketemplates(t_templateid(i)).template, 'EndTime', duration, 'SampleRate', sample_rate) ;
        if (i==1)
            spiketrainlength = length(targets(i).sampletargets) ;
        end
        % generate the differential and the second differential of both of these.
        ssamples_1d = simplediff(targets(i).sampletargets, dslength) ;
        ssamples_1dd = simplediff(ssamples_1d, dslength) ;
        % and normalise them
        targets(i).ssamples_1n = setlimits(targets(i).sampletargets, -0.5, 0.5) ;
        targets(i).ssamples_1dn = setlimits(ssamples_1d, -0.5, 0.5) ;
        targets(i).ssamples_1ddn = setlimits(ssamples_1dd, -0.5, 0.5) ;
        clear targettimes ;
    end
else
    % check that n_targets corresponds to the supplied target dataset.
    if (length(targets) ~= n_targets)
        error(['number of targets proposed differs from number of targets in the supplied dataset']) ;
    end
    spiketrainlength = length(targets(1).sampletargets) ;
end
% generate some jittered data
for nj = 1: n_jitter
    % generate jittered noisy spike trains
    % generate jittered spike trains from the appropriate target,
    % jitter_orig(nj)
    jstimes = probjitter(targets(jitter_orig(nj)).targettimes, jprobspikes(nj), jjstd(nj), duration );
    % then generate samples from the appropriate spike template,
    % j_templateid
    [jsamples temp] = gensampledspikes(jstimes, spiketemplates(j_templateid).template, 'EndTime', duration, 'SampleRate', sample_rate) ;
    if ~isequal(length(jsamples), spiketrainlength)
        disp(['Warning: jittered spike train length not same as target']) ;
    end
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
    % generate some sampled spike data
    [u_ssamples temp] = gensampledspikes(u_stimes, spiketemplates(u_templateid(nu)).template, 'EndTime', duration, 'SampleRate', sample_rate) ;
    if ~isequal(length(u_ssamples), spiketrainlength)
        disp(['Warning: uncorrelated spike train length not same as target']) ;
    end
    u_ssamples_d = simplediff(u_ssamples, dslength) ;
    u_ssamples_dd = simplediff(u_ssamples_d, dslength) ;
    % and normalise them
    u_samples_n(nu,:) = setlimits(u_ssamples, -0.5, 0.5) ;
    u_samples_dn(nu,:) = setlimits(u_ssamples_d, -0.5, 0.5) ;
    u_samples_ddn(nu,:) = setlimits(u_ssamples_dd, -0.5, 0.5) ;
    clear u_stimes ;
end

% save a little space
clear u_ssamples u_ssamples_d u_ssamples_dd 
% now mix these together
corrspikenoise = zeros([1 spiketrainlength] );
for i_corr = 1:n_jitter
    corrspikenoise = corrspikenoise + jmixcoeffs_orig(i_corr) * jsamples_n(i_corr, :) ;
    corrspikenoise = corrspikenoise + jmixcoeffs_d(i_corr) * jsamples_dn(i_corr, :) ;
    corrspikenoise = corrspikenoise + jmixcoeffs_dd(i_corr) * jsamples_ddn(i_corr, :) ;
end
% save some more space
clear jsamples_n jsamples_dn jsamples_ddn

uncorrspikenoise = zeros([1 spiketrainlength] );
for i_uncorr=1:n_uncorr
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_orig(i_uncorr) * u_samples_n(i_uncorr,:) ;
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_d(i_uncorr) * u_samples_dn(i_uncorr,:) ;
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_dd(i_uncorr) * u_samples_ddn(i_uncorr,:) ;
end

% clear some space
clear u_samples_n u_samples_dn u_samples_ddn

% add up the target (or original) signal
origsignal = zeros([1 spiketrainlength] );
for i_target = 1:n_targets
    origsignal = origsignal +  origmixcoeffs_orig(i_target) * targets(i_target).ssamples_1n +  ...
    origmixcoeffs_d(i_target) * targets(i_target).ssamples_1dn ...
    + origmixcoeffs_dd(i_target) * targets(i_target).ssamples_1ddn ; 
end

% add up the original, jittered and uncorrelated signals
noisysignal = origsignal + corrspikenoise + uncorrspikenoise ; 
% noisysignal = origsignal ;
% add some white noise
almostfinalsignal = awgn(noisysignal, noise_snr, 'measured') ;

% finally, scale the final signal 
finalsignal = setlimits(almostfinalsignal, finalminval, finalmaxval) ;

end

