function [finalsignal, targets, rngreturnstatus] = generatenoisysamples(varargin)
%
% This is version 1.12, released 27 July 2017
% 
% changes from 1.11:
% All calls to exist() now have 2 paramters (fixes issue when there's a
% directory locally with the same name as the variable
% All calls to interp1 now use pchip instead of cubic
%
% Has spatial effect, as 0.3, but this time effect
% different for orig, 1, and 2 deriv. Put the values out to separate
% functions
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
% Note that the spike times are the *starting times* of the spikesampletarge
% templates, so that the actual peak will be some little time later. This
% last set are included so that one can re-do data sets with the same
% primary spikes, but different noise levels etc. See varargin below.
%
% varargin: lots of them 
% N_Targets, followed by the number of target neurons
% ReuseTargets, followed by a structure which contains the actual times of
% spike starts (so that one can generate a number of sets of spikes with
% different noise levels etc., but the *same* underlying sets of primary
% spikes
% Duration, followed by the duration of the spike train (defaults to
% 0.05s)
% SampleRate, followed by the sample rate (defaults to 100000)
% DSLength, used in smoothing (defaults to 15 samples)
% ReuseRNGstate: re use the same random number generator state
%
% New version started 23 Oct 2005: have a number of target neurons, each with
% their own parameters, so that one can do spike sorting as well as spike
% detection. Further, allow the correlated spike trains to be correlated to
% one of the target neurons, and allow them to use whatever original spike
% train is required.
%
% new version started 29 Jan 2006: uses gensampledspikes_smoothed, and
% produces traces which take account of the spatial extent, and the resultant
% temporally extended signal.
%
% extended 
% LSS started 8 June 2005
% First release 14 June 2005.
% Added some error traps (to try to avoid silent crashes due to
% inappropriate parameter selection): 3 Nov 2005. LSS
% modifications: temporal extension: 29 Jan 2006 LSS
% Various fixes and improvements: 20 March 2006: LSS
% more of the same June 2006.
% minimal modification of file names, also made noise_snr in varargin
% LSS 29 2 2012 (leap!)
%
% 
% defaults:  may be overwritten in varargin section

n_targets = 2 ; % number of target neurons
n_jitter =  7; %number of jittered versions of the original spike train
n_uncorr = 15 ; % number of uncorrelated spike trains

duration = 0.5 ; % length of the signal: NB: approximate. actual length may be a little longer
sample_rate = 24000 ;
dslength = 12 ; % was 75 used in smoothing during diferentiation. Is a number of samples.
refractoryperiod = 0.001 ; % refractory period in seconds

% temporal/spatial extension parameters. % added 29 1 2006 LSS
t_delta_integrate = 0.000048 ; % delay between time integrate points 
% in seconds: 30 microseconds. Note that gensampledspikes_smoothed
% requires ths value in milliseconds.
n_delta_integrate = 25 ; % number of time integrate points
% location of the temporal integration values
temporaldatainfile = true; ;

% modified 29 2 2012: rhis way the data is in a subfolder example6,
% wherever the actual matlab file is loaded. If you want to move it, alter
% the line below.
% directory name for loading the temporal / spatial data from
% temporaldirectory = '/Users/lss/matlab_stuff/NoiseModellingv1point0/example6' ;
temporaldirectory = 'example6' ;

% noise level defaults to -100dB down (invisible!)
% in dB (for awgn function)
noise_snr = 100 ;

showsnr = 0 ; % default is not to show SNR
sametargetsizes = 0 ; % default is not to make target sizes all the same
t1wtssupplied = 0 ;
% overwrite earlier parameters if required
for i=1:(nargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'N_Targets'
                n_targets = varargin{i+1};
            case 'N_Jitter'
                n_jitter  = varargin{i+1};
            case 'N_Uncorr'
                n_uncorr = varargin{i+1};
            case 'ReuseTargets'
                targets = varargin{i+1};
            case 'Duration'
                duration = varargin{i+1};
            case 'SampleRate'
                sample_rate = varargin{i+1};
            case 'DSLength'
                dslength = varargin{i+1};
            case 'ReuseRNGState'
                rngstate = varargin{i+1};
            case 'RefractoryPeriod'
                refractoryperiod = varargin{i+1};
            case 'T_Delta_Integrate'
                t_delta_integrate = varargin{i+1};
            case 'N_Delta_Integrate'
                n_delta_integrate = varargin{i+1};
            case 'ReuseTargetTimes'
                oldtargettimes = varargin{i+1};
            case 'TemporalDataInFile'
                temporaldatainfile = varargin{i+1};
            case 'TemporalDirectory'
                temporaldirectory = varargin{i+1};
            case 'ShowSNR'
                showsnr = varargin{i+1};
            case 'Orig_Temporal_Supplied'
                orig_temporal = varargin{i+1};
            case 'SameTargetSizes'
                sametargetsizes = varargin{i+1};
            case 'Target1weights'
                t1wtssupplied = 1 ;
                t1wts = varargin{i+1};
            case 'NoiseSNR'
                noise_snr = varargin{i+1};
        end
    end
end

% templates for spike for target (primary) neuron. 
% File names no longer need to be same length (now using a cell array)
templatenames = {'spike_1n.csv',  
    'spike_2n.csv',
    'naundorf_vivo.csv'}; 
n_templates = length(templatenames) ; 

% spatial/temporal delay information: added 29 1 2006 LSS
% consists of a number of arrays: 
% if temporaldatainfile then the values will be read from a file (and this
% can be created using the PrepareParams Mac OSX application: otherwise the
% values are in the file themselves (see below)
% each array should be n_delta_integrate long, and there should be one for
% each neuron.
% spatiotemporal information for target neurons 
if temporaldatainfile
    if not(exist('orig_temporal', 'var'))
        orig_temporal = target_temporal_mac(n_targets, n_delta_integrate, temporaldirectory) ;
    end
    % spatiotemporal information for correlated noise neurons
    corr_temporal =  correlated_temporal_mac(n_jitter, n_delta_integrate, temporaldirectory ) ;
    % spatiotemporal informqtion for uncorrelated noise neurons
    n_uncorr_supplied = 15 ; % for now: if we alter uncorrelated_temporal, we may want to alter this too
    % note that dirstring is not actually used here.
    u_temporal = uncorrelated_temporal_mac(n_uncorr_supplied, n_delta_integrate, temporaldirectory) ;
else
    if not(exist('orig_temporal', 'var'))
        orig_temporal = target_temporal(n_targets, n_delta_integrate) ;
    end
    % spatiotemporal information for correlated noise neurons
    corr_temporal =  correlated_temporal(n_jitter, n_delta_integrate) ;
    % spatiotemporal informqtion for uncorrelated noise neurons
    n_uncorr_supplied = 15 ; % for now: if we alter uncorrelated_temporal, we may want to alter this too
    % note that dirstring is not actually used here.
    u_temporal = uncorrelated_temporal(n_uncorr_supplied, n_delta_integrate) ;
end


% target neuron information
% needs to be replicated n_targets times
% below is for n_targets = 0 to 4.

spikedist = ['P' ,% distribution to be used in generation:'P' for poission, 'G' for Gaussian
    'G',
    'P',
    'G', 
    'P',
    'G',
    'P',
    'G',
    'P'] ; 
% for poission distribution only
poissionnumber = [4 0 4 0 4 0 4 0 8] ; % poisson values: only used where there's a 'P' in spikedist
poissmeanISI = [0.01 0 0.07 0 0.08 0 0.07 0 0.015];


%for gaussian only
gaussmeanITD = [0 0.05  0 0.5 0 0.5 0 0.03 0]; % gaussian values: only used where there's a 'G' in spike dist
gaussstdev = [0 0.012 0 0.13 0 0.14 0 0.005 0 ] ;

t_templateid = [3 3 3 3 3 3 3 1 1] ; % template for this spike train
if (n_targets > 0 ) && (max(t_templateid) > n_templates) % can't use non-existant templates for spike generation
    error(['no of templates = ' num2str(n_templates) ', but highest template for target neuron = ', num2str(max(t_templateid)) ] ) ;
end
% for differential amounts of each of the signals
% mixing coefficients for the primary neural signal.
% morphing values
origmixcoeffs_orig = [0.075 0.05 1 1]; % plain (intracellular) signal strength
origmixcoeffs_d = [-0.35 0.33 1 1]; % differentiated signal strength
origmixcoeffs_dd = [0.04 0.1 1 1]; % twice differentiated signal strength
% origmixcoeffs_orig = [0.1 0.1 1 1]; % plain (intracellular) signal strength
% origmixcoeffs_d = [-0.1 -0.1 1 1]; % differentiated signal strength
% origmixcoeffs_dd = [0.05 0.05 1 1]; % twice differentiated signal strength
if t1wtssupplied > 0 % target 1 weights supplied (useful for morphing)
    origmixcoeffs_orig(1) = t1wts(1) ;
    origmixcoeffs_d(1) = t1wts(2) ;
    origmixcoeffs_dd(1) = t1wts(3) ;
end

% for jittered noise
[n_jitter_supplied temp] = size(corr_temporal) ; % number of jittered (correlated) spike trains supplied: loop around
% in varargin
jitter_orig = [1 1 1 1 1 1 1] ; % set all to 1's so it works with a single target neuron
% jitter_orig = [1 2 1 2 1 1 1] ;
% could be e.g. [1 2 3 3 4 4 1]  for use with 4 neurons; % which original signal should be jittered
if (n_jitter > 0 ) && (max(jitter_orig) > n_targets) % can't jitter targets that don't exist
    error(['no of targets = ' num2str(n_targets) ', but highest target neuron to jitter = ', num2str(max(jitter_orig)) ] ) ;
end
jitteroveralllevel = 2; %sqrt(2) ;
% jitteroverallsize_orig = 0.01 ; %overall size of the original jittered signal
% jitteroverallsize_d = 0.1 ; %overall size of the differentiated jittered signal
% jitteroverallsize_dd = 0.1 ; %overall size of the twice differentiated jittered signal
jitteroverallsize_orig = 0.008 * jitteroveralllevel; %overall size of the original jittered signal
jitteroverallsize_d = -0.08 * jitteroveralllevel; %overall size of the differentiated jittered signal
jitteroverallsize_dd = 0.02 * jitteroveralllevel; %overall size of the twice differentiated
% jittered signal
% The arrays below must ne at least n_jitter long
j_templateid = [1 2 3 1 2 3 1 ] ; % templates for these spikes 
if (n_jitter > 0 ) && (max(j_templateid) > n_templates) % can't jitter using templates that don't exist
    error(['no of templates = ' num2str(n_templates) ', but highest jitter neuron template = ', num2str(max(j_templateid)) ] ) ;
end
jitterdistribution = ['e' 'e' 'e' 'e' 'e' 'e' 'e'] ; % e for even, n for normal
jprobspikes = [0.6 0.8 0.7 0.5 0.3 0.5 0.7] ; % probability that an inpout spike results in an output spike
jjdelta = [0.01 0.008 0.01 0.012 0.011 0.009 0.01] ; % if jitterdistribution = 'n', 
% these are standard deviation of the jitter (mean = 0), and if
% jitterdistribution = e, these are the maximal delta times
%
% jittered noise mixing coefficients
jmixcoeffs_orig = [0.1 0.2 0.27 0.153 0.2 0.04 0.24 ] * jitteroverallsize_orig  ; % plain (intracellular signal strength
jmixcoeffs_d = [0.1 0.5 0.54 0.46 0.14 0.38 0.32] * jitteroverallsize_d ; % for 1st differential
jmixcoeffs_dd = [0.19 0.13 0.41 0.1 0.2 0.39 0.16] * jitteroverallsize_dd; % for 2nd differential

% for uncorrelated spike noise

ucoveralllevel = 24; %sqrt(2) ; 
% u_overallsize_orig = 0.004 ; %overall size of the original un-jittered signal
% u_overallsize_d = -0.04 ; %overall size of the differentiated un-jittered signal
% u_overallsize_dd = -0.02 ; %overall size of the twice differentiated un-jittered signal
u_overallsize_orig = 0.0016 * ucoveralllevel; %overall size of the original un-jittered signal
u_overallsize_d = -0.016 * ucoveralllevel; %overall size of the differentiated un-jittered signal
u_overallsize_dd = -0.008 * ucoveralllevel; %overall size of the twice differentiated un-jittered signal
% u_overallsize_orig = 0.0008 ; %overall size of the original un-jittered signal
% u_overallsize_d = -0.008 ; %overall size of the differentiated un-jittered signal
% u_overallsize_dd = -0.004 ; %overall size of the twice differentiated un-jittered signal

% All arrays below must be at least n_uncorr_supplied long
u_templateid = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 ]; % templates for these spikes 
if (n_uncorr > 0) && (max(u_templateid) > n_templates) % can't generate uncorrelated spikes using templates that don't exist
    error(['no of templates = ' num2str(n_templates) ', but highest uncorrelated neuron template = ', num2str(max(u_templateid)) ] ) ;
end
u_spikedist = ['P' 'G' 'G' 'G' 'P' 'P' 'G' 'G' 'G' 'P' 'P' 'G' 'G' 'G' 'P'] ; % distribution type for these spike trains
u_poissionnumber = [6 0 0 0 8 6 0 0 0 9 6 0 0 0 8] ; % number for use in generating Poisson distribution
u_poissmeanISI = [0.04 0 0 0  0.025 0.03 0 0 0  0.025 0.03 0 0 0  0.025] ; % Mean ISI for each poisson distribution
u_gaussmeanITD = [0 0.03 0.025 0.03 0 0 0.03 0.025 0.03 0 0 0.04 0.025 0.03 0] ; % Mean Interspike time difference (shoud be ISI!)
u_gaussstdev = [0 0.004 0.005 0.005 0 0 0.004 0.005 0.005 0  0 0.004 0.005 0.005 0 ] ; % STD thereof
% uncorrelated spikes mixing coefficients
umixcoeffs_orig = [0.133 0.022 0.139 0.012 0.122 0.33 0.22 0.16 0.12 0.26 0.34 0.29 0.39 0.17 0.20] * u_overallsize_orig ; % for intracellular signal
umixcoeffs_d = [0.337 0.45 0.59 0.21 0.84 0.33 0.44 0.59 0.21 0.54 0.33 0.41 0.53 0.26 0.64] * u_overallsize_d ; % for 1st differential
umixcoeffs_dd = [0.78 0.67 0.39 0.12 0.33 0.47 0.82 0.23 0.62 0.38 0.76 0.19 0.59 0.29 0.71] * u_overallsize_dd ; % for 2nd differential



% final output range
finalminval = -1.0 ;
finalmaxval = 1.0 ;

% end of parameters!


% read in the spike templates (used for targets and uncorrelated spikes)
for i = 1:n_templates
    spiketemplates(i).template = csvread(templatenames{i}) ;
    % adjust so that the times are norm alised to start at 0
    [tlength temp] = size(spiketemplates(i).template) ;
    tadj = ones([tlength 1]) * spiketemplates(i).template(1,1) ;
    spiketemplates(i).template(:,1) = spiketemplates(i).template(:,1) - tadj ;
    spiketemplates(i).templatename = templatenames(i,:) ;
end

% reset the random number generator if appropriate
if exist('rngstate', 'var') == 1
    rand('state',rngstate) ;
else
    rand('state',sum(100*clock)) ;
end

rngreturnstatus = rand('state') ; % in case we might want to use it again

% if targets has been supplied as varargin, then use them. Otherwise
% generate the target spike times. 
if ~exist('targets', 'var')
    % generate the target spike times
    for i = 1:n_targets
        if exist('oldtargettimes', 'var') == 1
            % re-use the target times
            if length(oldtargettimes) < i
                error(['re-using old target times, but there are more targets in this run']) ;
            end
            targettimes = oldtargettimes(i).targettimes ;
        else
            if spikedist(i) == 'P'
                targettimes = genspikespoisson(duration, poissionnumber(i), poissmeanISI(i), 'MinISI', refractoryperiod) ;
            else if spikedist(i) == 'G'
                    targettimes = genspikesgaussian(duration, gaussmeanITD(i), gaussstdev(i), 'MinISI', refractoryperiod) ;
                end
            end
        end
        targets(i).targettimes = targettimes ;

        % generate some sampled spike data from template t_templateid(i)
        % (not spatiotemporal: add the spatiotemporal effects in later.)
        [targets(i).sampletargets targets(i).actualpeaks] = ...
            gensampledspikes(targettimes, spiketemplates(t_templateid(i)).template, ...
            'EndTime', duration, 'SampleRate', sample_rate) ;

        
        % generate the differential and the second differential of both of these.
        ssamples_1d = simplediff(targets(i).sampletargets, dslength) ;
        ssamples_1dd = simplediff(ssamples_1d, dslength) ;
        % now perform the spatiotemporal transform on all three of these
        targets(i).stsamples = spatiotransform(targets(i).sampletargets,  sample_rate,  ...
                    'DelayDelta', t_delta_integrate, 'DelayNumber', n_delta_integrate, ...
                    'DelayCoefficients', squeeze(orig_temporal(i,1,:)) ) ;
        if (i==1)
            spiketrainlength = length(targets(i).stsamples) ;
        end
        
        st_ssamples_1d = spatiotransform(ssamples_1d,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(orig_temporal(i,2,:)) ) ;
        st_ssamples_1dd = spatiotransform(ssamples_1dd, sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(orig_temporal(i,3,:)) ) ;
        % and normalise them
        targets(i).ssamples_1n = setlimits(targets(i).stsamples, -0.5, 0.5) ;
        targets(i).ssamples_1dn = setlimits(st_ssamples_1d, -0.5, 0.5) ;
        targets(i).ssamples_1ddn = setlimits(st_ssamples_1dd, -0.5, 0.5) ;
        targets(i).final = origmixcoeffs_orig(i) * targets(i).ssamples_1n + ...
            origmixcoeffs_d(i) * targets(i).ssamples_1dn + origmixcoeffs_dd(i) * targets(i).ssamples_1ddn ;
        if sametargetsizes > 0
            targets(i).final = setlimits(targets(i).final, -0.1 , 0.1) ;
        end
        
    end
    clear targettimes ssamples_1d ssamples_1dd st_ssamples_1d st_ssamples_1dd;
else
    % check that n_targets corresponds to the supplied target dataset.
    if (length(targets) ~= n_targets)
        error(['number of targets proposed differs from number of targets in the supplied dataset']) ;
    end
    spiketrainlength = length(targets(1).stsamples) ;
end
% generate some jittered data
for nj_1 = 1: n_jitter
    % generate jittered noisy spike trains
    nj = mod(nj_1, n_jitter_supplied) ; % go round and round if required
    if (nj == 0) 
        nj = n_jitter_supplied ; 
    end 
    % generate jittered spike trains from the appropriate target,
    % note that these jittered times may be evenly or normally distributed
    % jitter_orig(nj)
    if (jitterdistribution(nj) == 'n')
        jstimes = probjitter(targets(jitter_orig(nj)).targettimes, jprobspikes(nj), jjdelta(nj), duration );
    else if (jitterdistribution(nj) == 'e')
            jstimes = probjitter_even(targets(jitter_orig(nj)).targettimes,...
                jprobspikes(nj), jjdelta(nj), duration );
        end
    end
    % then generate samples from the appropriate spike template,
    % j_templateid
    [jsamples temp] = gensampledspikes(jstimes, spiketemplates(j_templateid(nj)).template, ...
        'EndTime', duration, 'SampleRate', sample_rate) ;
%         'DelayDelta', t_delta_integrate * 1000, 'DelayNumber', n_delta_integrate, ...
%         'DelayCoefficients', corr_temporal(nj,:)) ;
    
    % and find its derivatives
    jsamples_d = simplediff(jsamples, dslength) ;
    jsamples_dd = simplediff(jsamples_d, dslength) ;
    

    st_jsamples = spatiotransform(jsamples,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(corr_temporal(nj,1,:)) ) ;      
    if (~exist('spiketrainlength', 'var') & (nj==1))
            spiketrainlength = length(st_jsamples) ;
    end
    if ~isequal(length(st_jsamples), spiketrainlength)
        disp(['Warning: jittered spike train length not same as target']) ;
    end
    st_jsamples_1d = spatiotransform(jsamples_d,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(corr_temporal(nj,2,:)) ) ;  
    st_jsamples_2d = spatiotransform(jsamples_dd,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(corr_temporal(nj,3,:)) ) ;  
    % and normalise the derivatives
    jsamples_n(nj_1, :) = setlimits(st_jsamples, -0.5, 0.5) ;
    jsamples_dn(nj_1, :) = setlimits(st_jsamples_1d, -0.5, 0.5) ;
    jsamples_ddn(nj_1, :) = setlimits(st_jsamples_2d, -0.5, 0.5) ;

end
clear jstimes jsamples jsamples_d jsamples_dd st_jsamples st_jsamples_1d st_jsamples_2d ;


% generate some uncorrelated spikes too
% modded LSS 30 Jan 2006 to allow many more uncorrelated spike trains by 
% looping around the data supplied
% modded LSS 6 6 06 to accumulate rather than store and add later
uncorrspikenoise = zeros([1 spiketrainlength] );
for nu_1=1:n_uncorr 
    % set nu so that nu_1 mod n_uncorr_supplied = 0 turns into
    % n_uncorr_supplied
    nu = mod(nu_1, n_uncorr_supplied) ;
    if (nu == 0) 
        nu = n_uncorr_supplied ; 
    end   

    if u_spikedist(nu) == 'P'
        u_stimes = genspikespoisson(duration, u_poissionnumber(nu), u_poissmeanISI(nu), 'MinISI', refractoryperiod) ;
    else if u_spikedist(nu) == 'G'
            u_stimes = genspikesgaussian(duration, u_gaussmeanITD(nu), u_gaussstdev(nu), 'MinISI', refractoryperiod) ;
        end
    end
    % generate some sampled spike data
    if isempty(u_stimes)
        if ~exist('spiketrainlength', 'var')
            spiketrainlength = duration * sample_rate ; % defensive
        end
        u_ssamples = zeros([1 spiketrainlength]) ;  
    else
        [u_ssamples temp] = gensampledspikes(u_stimes, spiketemplates(u_templateid(nu)).template,...
            'EndTime', duration, 'SampleRate', sample_rate) ;
    end
    % generate the derivatives
    u_ssamples_d = simplediff(u_ssamples, dslength) ;
    u_ssamples_dd = simplediff(u_ssamples_d, dslength) ;
    
    % perform spatiotemporal filtering
    st_usamples = spatiotransform(u_ssamples,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(u_temporal(nu,1,:)) ) ;      
    if (~exist('spiketrainlength', 'var') & (nu_1==1))
            spiketrainlength = length(st_usamples) ;
    end
    if ~isequal(length(st_usamples), spiketrainlength)
        disp(['Warning: uncorrelated spike train length not same as target']) ;
    end
    st_usamples_1d = spatiotransform(u_ssamples_d,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(u_temporal(nu,2,:)) ) ;  
    st_usamples_2d = spatiotransform(u_ssamples_dd,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(u_temporal(nu,3,:)) ) ;  
        


    % and normalise them
    u_samples_n = setlimits(st_usamples, -0.5, 0.5) ;
    u_samples_dn = setlimits(st_usamples_1d, -0.5, 0.5) ;
    u_samples_ddn = setlimits(st_usamples_2d, -0.5, 0.5) ;
    % accumulate the values
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_orig(nu) * u_samples_n ;
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_d(nu) * u_samples_dn ;
    uncorrspikenoise = uncorrspikenoise + umixcoeffs_dd(nu) * u_samples_ddn ;
end

% save a little space
clear u_ssamples u_ssamples_d u_ssamples_dd  st_usamples st_usamples_1d st_usamples_2d ;
clear u_samples_n u_samples_dn u_samples_ddn
% now mix these together
corrspikenoise = zeros([1 spiketrainlength] );
for i_corr_1 = 1:n_jitter
    i_corr = mod(i_corr_1, n_jitter_supplied) ; % go round and round if required
    if (i_corr == 0) 
        i_corr = n_jitter_supplied ; 
    end 
    corrspikenoise = corrspikenoise + jmixcoeffs_orig(i_corr) * jsamples_n(i_corr_1, :) ;
    corrspikenoise = corrspikenoise + jmixcoeffs_d(i_corr) * jsamples_dn(i_corr_1, :) ;
    corrspikenoise = corrspikenoise + jmixcoeffs_dd(i_corr) * jsamples_ddn(i_corr_1, :) ;
end
% save some more space
clear jsamples_n jsamples_dn jsamples_ddn

% 
% for i_uncorr=1:n_uncorr
%     nu = mod(i_uncorr, n_uncorr_supplied) ;
%     if (nu == 0) 
%         nu = n_uncorr_supplied ; 
%     end   
%     uncorrspikenoise = uncorrspikenoise + umixcoeffs_orig(nu) * u_samples_n(i_uncorr,:) ;
%     uncorrspikenoise = uncorrspikenoise + umixcoeffs_d(nu) * u_samples_dn(i_uncorr,:) ;
%     uncorrspikenoise = uncorrspikenoise + umixcoeffs_dd(nu) * u_samples_ddn(i_uncorr,:) ;
% end

% clear some space


% add up the target (or original) signal
origsignal = zeros([1 spiketrainlength] );
for i_target = 1:n_targets
%     origsignal = origsignal +  origmixcoeffs_orig(i_target) * targets(i_target).ssamples_1n +  ...
%     origmixcoeffs_d(i_target) * targets(i_target).ssamples_1dn ...
%     + origmixcoeffs_dd(i_target) * targets(i_target).ssamples_1ddn ; 
    origsignal = origsignal + targets(i_target).final ;
end

if (showsnr > 0)
    origpeak = max(abs(origsignal - mean(origsignal))) ;
    corrpeak = max(abs(corrspikenoise - mean(corrspikenoise))) ;
    uncorrpeak = max(abs(uncorrspikenoise - mean(uncorrspikenoise))) ;
    disp(['Orig peak = ' num2str(origpeak) ' Corr peak = ' num2str(corrpeak) ' Uncorr peak = ' num2str(uncorrpeak)]) ;
    origmedian = median(abs(origsignal - mean(origsignal))) ;
    corrmedian = median(abs(corrspikenoise - mean(corrspikenoise))) ;
    uncorrmedian = median(abs(uncorrspikenoise - mean(uncorrspikenoise))) ;
    disp(['Orig median = ' num2str(origmedian) ' Corr median = ' num2str(corrmedian) ' Uncorr median = ' num2str(uncorrmedian)]) ;
    origmean = mean(abs(origsignal - mean(origsignal))) ;
    corrmean = mean(abs(corrspikenoise - mean(corrspikenoise))) ;
    uncorrmean = mean(abs(uncorrspikenoise - mean(uncorrspikenoise))) ;
    disp(['Orig mean = ' num2str(origmean) ' Corr mean = ' num2str(corrmean) ' Uncorr mean = ' num2str(uncorrmean)]) ;
    % origstd = std((origsignal - mean(origsignal))) ;
    corrstd = std(corrspikenoise) ;
    uncorrstd = std(uncorrspikenoise) ;
    totalstd = std(corrspikenoise + uncorrspikenoise) ;
    disp([' Corr std = ' num2str(corrstd) ' Uncorr std = ' num2str(uncorrstd) 'Total std = ' num2str(totalstd)]) ;
end

% add up the original, jittered and uncorrelated signals
noisysignal = origsignal + corrspikenoise + uncorrspikenoise ; 
% noisysignal = origsignal ;
% add some white noise
almostfinalsignal = awgn(noisysignal, noise_snr, 'measured') ;

% finally, scale the final signal 
finalsignal = setlimits(almostfinalsignal, finalminval, finalmaxval) ;

if ~exist('targets', 'var')
    targets = [] ; % if there's only noise being generated, return a null satructure
end
end

