function  testthrl(ntg)
if ~exist('ntg', 'var') ntg=3; end
[setsignal targets] = gensampl('N_Targets', ntg,  'N_Jitter', 0,  'N_Uncorr', 2, 'Duration', 5, 'SampleRate', 24e3);
noisysignal= setsnr(setsignal, -10);
[pmid event c_spike peaktrain] = estcob(noisysignal(:,1), targets, 24e3, 1, 1, .5e-3);
j=sort(peaktrain); 
if j(end)-j(end-10)>=.8 
    peaktrain(peaktrain>=j(end-10))=j(end-10); 
end
peaktrain=peaktrain/max(abs(peaktrain));
num=[]; thrl99 = max(peaktrain)-max(peaktrain)*.01; stp=-max(peaktrain)*.01;
for thrl=thrl99:stp:0
    spikes= find(peaktrain>thrl);
    sptrain=[];  sptrain(spikes)=1;
    id=find(sptrain>0);
    difid=diff(id);
    l=0; r=id(1);
    for n=1:length(difid)
        if difid(n)<12
            l=l+difid(n);
            r(end)=floor(mean([id(n)+l, r(end)]));
        else
            r(n)=id(n+1); l=0;
        end
    end
    res=r(r>0);

    num=[num; [thrl length(res)]];
end
figure; bar(num(1:99,1), num(1:99,2));
disp([pmid])
if ntg==2
    disp([length(targets(1,1).targettimes), length(targets(1,2).targettimes)])
    disp([length(targets(1,1).targettimes)+ length(targets(1,2).targettimes)])
elseif ntg==3
    length ([targets(1,1).targettimes, targets(1,2).targettimes, targets(1,3).targettimes])
    disp([length(targets(1,1).targettimes), length(targets(1,2).targettimes), length(targets(1,3).targettimes)])
end

 
 
%**************************************************************************
%                                Function
%**************************************************************************
function [setsignal targets rngreturnstatus] = gensampl(varargin)
% Modified function of generatnoisysamples by LSS
%function [finalsignal targets rngreturnstatus] = generatenoisysamples(varargin)
%
% Default ==============================================================
% function [setsignal targets rngreturnstatus] = gensamples('N_Targets', 2, ...
%     'N_Jitter', 7,  'N_Uncorr', 15, 'Duration', 5, 'SampleRate', 24e3,  'DSLength', 12,  ...
%     'RefractoryPeriod', 0.001, 'T_Delta_Integrate', 0.000048, 'N_Delta_Integrate', 25, ...
%     'TemporalDataInFile', true,  'ShowSNR', 0, 'SameTargetSizes', 0)
%
%  'ReuseTargets', 'ReuseRNGState','TemporalDirectory','ReuseTargetTimes',
%  'Orig_Temporal_Supplied',  'Target1weights',
%=====================================================================
% This is version 1.1, released 20 June 2006.
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
temporaldatainfile = true; 

% directory name for loading the temporal / spatial data from
% temporaldirectory = 'C:\Documents and Settings\ssh\My Documents\Mtetwa_Leslie\NoiseModellingv1point1\NoiseModellingv1point1\example6'; %'H:\NoiseModellingv1point1\NoiseModellingv1point1\example6' ;
temporaldirectory = '/Users/lss/matlab_stuff/NoiseModellingv1point0/example6' ;
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
    if not(exist('orig_temporal'))
        orig_temporal = target_temporal_mac(n_targets, n_delta_integrate, temporaldirectory) ;
    end
    % spatiotemporal information for correlated noise neurons
    corr_temporal =  correlated_temporal_mac(n_jitter, n_delta_integrate, temporaldirectory ) ;
    % spatiotemporal informqtion for uncorrelated noise neurons
    n_uncorr_supplied = 15 ; % for now: if we alter uncorrelated_temporal, we may want to alter this too
    % note that dirstring is not actually used here.
    u_temporal = uncorrelated_temporal_mac(n_uncorr_supplied, n_delta_integrate, temporaldirectory) ;
else
    if not(exist('orig_temporal'))
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
poissmeanISI = [0.05 0 0.07 0 0.08 0 0.07 0 0.015];
%poissmeanISI = [0.013 0 0.07 0 0.08 0 0.07 0 0.015];
%poissmeanISI = [0.01 0 0.07 0 0.08 0 0.07 0 0.015];


%for gaussian only
%gaussmeanITD = [0 0.05  0 0.5 0 0.5 0 0.03 0]; % gaussian values: only used where there's a 'G' in spike dist
gaussmeanITD = [0 0.025  0 0.5 0 0.5 0 0.03 0]; % gaussian values: only used where there's a 'G' in spike dist
gaussstdev = [0 0.012 0 0.13 0 0.14 0 0.005 0 ] ;

t_templateid = [3 3 3 3 3 3 3 1 1] ; % template for this spike train
if (n_targets > 0 ) && (max(t_templateid) > n_templates) % can't use non-existant templates for spike generation
    error(['no of templates = ' num2str(n_templates) ', but highest template for target neuron = ', num2str(max(t_templateid)) ] ) ;
end
% for differential amounts of each of the signals
% mixing coefficients for the primary neural signal.
% morphing values
origmixcoeffs_orig = [0.075 0.05 .1 1]; % plain (intracellular) signal strength
origmixcoeffs_d = [-0.35 0.33 .1 1]; % differentiated signal strength
origmixcoeffs_dd = [0.04 0.1 0.11 1]; % twice differentiated signal strength
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
%jitter_orig = [1 2 1 2 1 1 1] ;
% could be e.g. [1 2 3 3 4 4 1]  for use with 4 neurons; % which original signal should be jittered
if (n_jitter > 0 ) && (max(jitter_orig) > n_targets) % can't jitter targets that don't exist
    error(['no of targets = ' num2str(n_targets) ', but highest target neuron to jitter = ', num2str(max(jitter_orig)) ] ) ;
end
jitteroveralllevel = 4; % sqrt(2) ;
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

ucoveralllevel = 4; %* sqrt(2) ; 
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

% noise level
% in dB (for awgn function)
noise_snr = 100 ;

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
if exist('rngstate') == 1
    rand('state',rngstate) ;
else
    rand('state',sum(100*clock)) ;
end

rngreturnstatus = rand('state') ; % in case we might want to use it again

% if targets has been supplied as varargin, then use them. Otherwise


%===============================================================================
% generate the target spike times. 
if ~exist('targets')
    % generate the target spike times
    for i = 1:n_targets
        if exist('oldtargettimes') == 1
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

%===============================================================================
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
    [jsamples temp] = gensampledspikes(jstimes, spiketemplates(j_templateid(i)).template, ...
        'EndTime', duration, 'SampleRate', sample_rate) ;
%         'DelayDelta', t_delta_integrate * 1000, 'DelayNumber', n_delta_integrate, ...
%         'DelayCoefficients', corr_temporal(nj,:)) ;
    
    % and find its derivatives
    jsamples_d = simplediff(jsamples, dslength) ;
    jsamples_dd = simplediff(jsamples_d, dslength) ;
    

    st_jsamples = spatiotransform(jsamples,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(corr_temporal(i,1,:)) ) ;      
    if (~exist('spiketrainlength') & (nj==1))
            spiketrainlength = length(st_jsamples) ;
    end
    if ~isequal(length(st_jsamples), spiketrainlength)
        disp(['Warning: jittered spike train length not same as target']) ;
    end
    st_jsamples_1d = spatiotransform(jsamples_d,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(corr_temporal(i,2,:)) ) ;  
    st_jsamples_2d = spatiotransform(jsamples_dd,  sample_rate, ...
            'DelayDelta', t_delta_integrate , 'DelayNumber', n_delta_integrate, ...
            'DelayCoefficients', squeeze(corr_temporal(i,3,:)) ) ;  
    % and normalise the derivatives
    jsamples_n(nj_1, :) = setlimits(st_jsamples, -0.5, 0.5) ;
    jsamples_dn(nj_1, :) = setlimits(st_jsamples_1d, -0.5, 0.5) ;
    jsamples_ddn(nj_1, :) = setlimits(st_jsamples_2d, -0.5, 0.5) ;

end
clear jstimes jsamples jsamples_d jsamples_dd st_jsamples st_jsamples_1d st_jsamples_2d ;


%===============================================================================
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
    if length(u_stimes) == 0
        if ~exist('spiketrainlength')
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
    if (~exist('spiketrainlength') & (nu_1==1))
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


%===============================================================================
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

% if (showsnr > 0)
%     origpeak = max(abs(origsignal - mean(origsignal))) ;
%     corrpeak = max(abs(corrspikenoise - mean(corrspikenoise))) ;
%     uncorrpeak = max(abs(uncorrspikenoise - mean(uncorrspikenoise))) ;
%     disp(['Orig peak = ' num2str(origpeak) ' Corr peak = ' num2str(corrpeak) ' Uncorr peak = ' num2str(uncorrpeak)]) ;
%     origmedian = median(abs(origsignal - mean(origsignal))) ;
%     corrmedian = median(abs(corrspikenoise - mean(corrspikenoise))) ;
%     uncorrmedian = median(abs(uncorrspikenoise - mean(uncorrspikenoise))) ;
%     disp(['Orig median = ' num2str(origmedian) ' Corr median = ' num2str(corrmedian) ' Uncorr median = ' num2str(uncorrmedian)]) ;
%     origmean = mean(abs(origsignal - mean(origsignal))) ;
%     corrmean = mean(abs(corrspikenoise - mean(corrspikenoise))) ;
%     uncorrmean = mean(abs(uncorrspikenoise - mean(uncorrspikenoise))) ;
%     disp(['Orig mean = ' num2str(origmean) ' Corr mean = ' num2str(corrmean) ' Uncorr mean = ' num2str(uncorrmean)]) ;
%     % origstd = std((origsignal - mean(origsignal))) ;
%     corrstd = std(corrspikenoise) ;
%     uncorrstd = std(uncorrspikenoise) ;
%     totalstd = std(corrspikenoise + uncorrspikenoise) ;
%     disp([' Corr std = ' num2str(corrstd) ' Uncorr std = ' num2str(uncorrstd) 'Total std = ' num2str(totalstd)]) ;
% end

% add up the original, jittered and uncorrelated signals
noisysignal = origsignal + corrspikenoise + uncorrspikenoise ; 
% noisysignal = origsignal ;
% add some white noise
% almostfinalsignal = awgn(noisysignal, noise_snr, 'measured') ;

% finally, scale the final signal 
% finalsignal = setlimits(almostfinalsignal, finalminval, finalmaxval) ;
 finalsignal = setlimits(noisysignal, finalminval, finalmaxval) ;
 
if ~exist('targets')
    targets = [] ; % if there's only noise being generated, return a null satructure
end
setsignal=[noisysignal' origsignal' corrspikenoise' uncorrspikenoise'];

%**************************************************************************
%                                Function
%**************************************************************************
function noisysignal= setsnr(setsignal, desiresnr)

if size(setsignal,2)~=4 
    error('wrong format of setsignal');
end

if ~exist('desiresnr') desiresnr=10; end
% if desiresnr>100 | desiresnr<-5 
%     disp(['Desired SNR is beyond the limit [20 ~ -5']);
%     disp(['Running program with the desiresnr =10db']);
%     desiresnr=10;
% end

% if ~exist('setsignal')
%     [setsignal targets rngreturnstatus] = makesamples('N_Targets', 2, 'Duration', 5);
% end

%[setsignal targets rngreturnstatus] = makesamples('N_Targets', 2, 'Duration', 5);
origsignal=setsignal(:,2);
corrspikenoise=setsignal(:,3);
uncorrspikenoise=setsignal(:,4);
fact=1;

%aveorigsignal=sqrt(mean(origsignal.^2));
aveorigsignal=abs(max(origsignal)-min(origsignal));
noise=corrspikenoise+uncorrspikenoise;
%avenoise=sqrt(mean(noise.^2));
avenoise=abs(max(noise)-min(noise));


thissnr= 20*log10(aveorigsignal/avenoise);
if thissnr>(desiresnr+.5) | thissnr<(desiresnr-.5)
    fact=10^((thissnr-desiresnr)/20);
    avenoise=avenoise*fact;
    thissnr= 20*log10(aveorigsignal/avenoise);
end
setsignal(:, [3, 4])= fact*[setsignal(:,3) setsignal(:,4)];
setsignal(:, 1)= setsignal(:,2)+setsignal(:,3)+setsignal(:,4);

noisysignal=setsignal; 
%plot(noisysignal(:, 1));
return


%**************************************************************************
%                                Function
%**************************************************************************
function [pmid event spike peaktrain] = estcob(dataset, gt, fs, misspenalty, insertpenalty, mintime_m)
%Find least penalty & best associates by CoB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       CoB BASED DETECTION METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmid(4,1) = 0; fl=0;
best_cob = 100000 ;

% Number of Actual Peaks
groundtruth = [] ;  for i = 1:length(gt)    groundtruth = [groundtruth  gt(i).actualpeaks] ; end
groundtruth = sort(groundtruth) ; ngt=length(groundtruth);
kact= floor(groundtruth*fs);
event(length(dataset),1)=0; event(kact)=1;  event=event(1:length(dataset));

[peaktrain]=spbycob(dataset, 256); peaktrain=peaktrain/max(abs(peaktrain));

% **********find best result************
tflag= true; thrl = .0;
while (tflag && (thrl < 3) && (thrl >=.0))
    thrl=thrl+0.01;
    %=========   thresholding ===========
    spikes= find(peaktrain>thrl);
    if ~isempty(spikes); c_spike = zeros(length(dataset), 1); c_spike (spikes)= 1; end;
    kest1=find(c_spike(1:end)>0); id=find(diff(kest1)< (mintime_m*fs));
    %================Estimate Spikes, Missing & Insertion===========
    
    %c_spike =findimpulse(dataset, 256, thrl); kest1=[find(c_spike(124:end)>0)]; id=find(diff(kest1)< (mintime_m*fs));
    %c_spike =neuron_impulse(dataset, 256, thrl);kest1=[find(c_spike(124:end)>0)]; id=find(diff(kest1)< (mintime_m*fs));
    % Estimated Position of Peaks
    

    
    if ~isempty(id)
        temp=(1:id(1)); for i=2:length(id) temp=[temp 2+id(i-1):id(i)]; end
        temp=[temp 2+id(end): length(kest1)]; kest=kest1(temp);
    else kest=kest1;
    end

    % ***************for any mintime_m *********************

    %Estimate the location of the spikes in Seconds
    pos_len = length(kest);
    T_pos = zeros(1, pos_len);
    for i=1:pos_len
        T_pos(i) = kest(i) /fs;
    end

    % find number of missing spikes
    [spike_mat d_cob_miss] = matchspikes_sd(gt, kest,'MinTime', mintime_m, 'SampleRate', (fs)) ;
    %find inserted spikes
    d_cob_ins = findinsertions_sd2(T_pos, 'Matched', spike_mat, 'srate', (fs)) ;

    %====================================================
    % calculate penalty
    penalty = misspenalty * d_cob_miss + insertpenalty * length(d_cob_ins) ;

    %*******************************************************


    % ***************for mintime_m 1.25e-4*********************
    %
    %     % Find number of missing spike  at signal.
    %     match_sp=sort([find(ismember(kact, kest))  find(ismember(kact, kest+1))  find(ismember(kact, kest+2)) ...
    %         find(ismember(kact, kest-1))  find(ismember(kact, kest-2)) ]');
    %     d_cob_miss=kact'; d_cob_miss(match_sp)=[];
    %
    %     % Find number of inserted spike at Estimated spike Matrix -- kest.
    %     match_sp=sort([find(ismember(kest, kact))  find(ismember(kest, kact+1))  find(ismember(kest, kact+2)) ...
    %         find(ismember(kest, kact-1))  find(ismember(kest, kact-2)) ]');
    %     d_cob_ins = kest'; d_cob_ins(match_sp)=[];
    %
    %====================================================
    % calculate penalty
    % penalty = misspenalty * length(d_cob_miss )+ insertpenalty * length(d_cob_ins) ;

    %*******************************************************

    if penalty < best_cob
        best_cob = penalty ;  % store if best so far
        best_miss = d_cob_miss;
        best_insert =length(d_cob_ins);
        best_thrl=thrl;
        fl=0;    spike=c_spike; spike(:)=0; spike(kest)=1;
    else
        fl=fl+1;
        if (best_cob==0  || fl==5); tflag=false; end
    end

end
%pmid(1) = 100*best_cob/ngt;
pmid(1) = 100*(ngt-best_miss)/ngt;
pmid(2) = 100*best_miss/ngt;
pmid(3) = 100*best_insert/ngt;
pmid(4) = best_thrl;



%**************************************************************************
%                                Function
%**************************************************************************

function [peaktrain]=spbycob(signal, nfft)
%estimate spike event train by CoB (also spikes with noisy environment)

%  spikes = spbycob(signal, nfft, thrl, display)
%  find the spike index in data using
%
% Example: 
%  spikes = sp_cob(signal, 256, 1.2, true)
%
%   peaktrain -     Output >> a logical peak indication  of spike + noise event
%   signal - -      Signal (spike train to be analyzed)
%   nfft    -        number of fourier point

% ****************  EXAMPLE   ****************
%sp_train=[];
%for i=1: length(data)/25e3
%    signal=data((i-1)*25e3+[1:25e3],2); 
%    [train]=spbycob(signal, 128,  0, 0); 
%    sp_train=[sp_train train(:)];
%end;
%sp_train=sp_train(:);
% ****************EXAMPLE  END ****************



signal=signal-mean(signal);   signal=signal./max(abs(signal)); signal =signal(:);
lt=length(signal);
% ------------------ Inverse Filter -----------------------------
[iht]=ifilt(signal,nfft);

%----------------------Noisy Impulse find -----------------------
signal=[signal(1:end-1)+signal(2:end); 0]*.5;
ft=conv(signal,iht); 
ft=ft(:); ft=[zeros(nfft/2+6,1); ft((nfft+1):(end-nfft))]; 

%================Denoising =====================
xtra0=length(ft)-nfft*(floor(length(ft)/nfft)); 
if xtra0>0; ft=[ft; zeros(2^nextpow2(xtra0)-xtra0,1)]; end;    % Zero padding to use wevelet function

%[s] = wden(ft,'rigrsure','s','mln',n,'bior1.5');

temp=ft; [a,d]=swt(ft, 3, 'coif1');  
if abs(skewness(temp))<1
    for n=1:3
        s=iswt(a(n,:), zeros(1, length(d)), 'coif1')';
        if abs(skewness(s))> abs(skewness(temp))  
            temp=s;
            if abs(skewness(s))>=1; break; end
        end
    end
end

if abs(skewness(temp))<1
    for n=1:3
        s=iswt( zeros(1, length(d)), d(n,:), 'coif1')';
        if abs(skewness(s))> abs(skewness(temp))    
            temp=s; %disp('used D');
            if abs(skewness(s))>=1; break; end
        end
    end
end
ft=temp; 

%----------------------Noise Suppressing -------------
f(lt,1)=0;  f(1:length(ft))=ft;
t1=f.*(f>=0);                    %t1 matches at s(nfft/2+7:end -nfft/2-6) for WT &   s(nfft/2:end -nfft/2)
t2=t1.*t1;
t3=t2.*t1;g=sort(t3,'descend'); t=(10*t3/min(g(1:10)));


abovethrl= find(t>0); h=[abovethrl(1);  abovethrl(2:end).*(diff(abovethrl)~=1)]; pulses=h(h>0);
sp_id=[]; for i=2:length(pulses); [val ind]=max(t(pulses(i-1):pulses(i))); sp_id(i-1)=pulses(i-1)+ind-1; end
sp_id=[sp_id, pulses(end)];
peaktrain= zeros(lt, 1); peaktrain(sp_id)= t(sp_id);
return


%================== FUNCTION Inverse Filter Estimation ==================

function [ht1, flag]=ifilt(x,nfft)

flag=0;
[bx,bc]=bg(x,nfft);
if sum(real(bx(:)))<0
    [bx,bc]=bg(-x,nfft); flag=1;
end

lb=log(bx); b=ifft(lb);
hx=([b(1,1:end)]); cbhx=hx(2:end)-real(b(1,1)); cbhx=[hx]; cbhx=cbhx-mean(cbhx); chx=exp(cbhx);


% ----------- Fourier Phase Estimation ---------------------

N=nfft/2;

if rem(N,2)==0
    Nby2=N+1;
else
    Nby2=N;
end
r=1;
Nby4=fix(Nby2/2)+1;
R=2*sum(1:Nby4-1)+Nby4-Nby2+1;

psi=zeros(R,1);
amat=zeros(R,Nby2);
for k=2:Nby4
    for l=2:k
        r=r+1;
        psi(r,1)=angle(bx(l,k));
        amat(r,k)=amat(r,k)+1;
        amat(r,l)=amat(r,l)+1;
        amat(r,l+k-1)=amat(r,l+k-1)-1;
    end
end
for k=Nby4+1:Nby2-1
    for l=2:(Nby2+1-k)
        r=r+1;
        psi(r,1)=angle(bx(l,k));
        amat(r,k)=amat(r,k)+1;
        amat(r,l)=amat(r,l)+1;
        amat(r,l+k-1)=amat(r,l+k-1)-1;
    end
end

amat(1,1)=1; psi(1,1)=angle(bx(1,1));

%-------------------- phase unwrapping -----------------------
PHI=[angle(chx(1:end)) 0];
kwrap = fix( (amat*PHI' - psi(1:end)) /(2*pi));
phi = amat(:,1:Nby2-1) \ (psi(1:end) + 2*pi*kwrap);

%--------------- Magnitude & Phase for IFFT ------------------
mag=abs(chx);
mag=[mag]';  mag = [mag(1:Nby2-1); flipud(mag(1:Nby2-2))];
phz = [phi(1:Nby2-1); -flipud(phi(1:Nby2-2))];
chx =( 1./mag) .* exp(-sqrt(-1)*phz);


ht1=ifft(chx);
%ht1=(ifftshift(ht1));
ht1=real(ifftshift(ht1));
if flag==1  ht1=-ht1;  end

return


%=================== FUNCTION Bispectrum Estimation ===================

function [Bx,Bc, Px]=bg(data, nfft)


if isempty(data) error('No data to analyze'); end

[Row, Col]= size(data); if (min(Row, Col)~=1) data=data(:); end

len=length(data);
if (floor(len/nfft)==0) error('short data length'); end

%------------------------------- Data Segmentation in M Realization and K Samples each =======

samp=nfft;                % Samples per realization
realz=fix(length(data)/samp);                  % Number of realization
data=data(1:samp*realz);
data=reshape(data,samp,realz);


%------------------------------- Window setting -------------------------------

wind=hann(samp); wind=wind(:);

%------------------------------- Spectrum -------------------------------

if (rem(nfft,2)~=0)
    mrow =fix(nfft/2)+1;
else
    mrow=nfft/2;
end
ncol=mrow;

Px = zeros(nfft,1);
Bx = zeros(mrow,ncol);

mask = hankel([1:mrow],[mrow:mrow+ncol-1] );   % the hankel mask (faster)

Nsamp = [1:samp]';
for Nrealz = 1:realz
    xseg = data(Nsamp)- mean(data(Nsamp)); % Subtract mean value from each record
    xseg   = xseg.*wind;                   % Passing data through window
    Xf     = fft(xseg, nfft)/samp;
    CXf    = conj(Xf);

    Bx  = Bx + (Xf(1:mrow) * Xf(1:ncol).').* CXf(mask);
    Px = Px + Xf.*CXf;

    Nsamp = Nsamp + samp;
end

Px = Px/(realz);        % averaging
Bx = Bx/(realz);        % Bispectrum averaging

%------------------------------- find out the Bicoherence -------------------------------

nm=(Px(1:mrow)*Px(1:ncol).').*Px(mask);
Bc=Bx.^2./(nm);

return
