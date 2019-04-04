function finalnewtimes = probjitter_even(oldtimes, probspike, jitterlength, duration)
%
% produces a new list of times from the old list. Each output spike is
% produced from each input spike with probability probspike, at a time
% randomly chosen between the t_orig - jitterlength and t_orig +
% jitterlength
%
% no additional spikes cuirrently produced.
%
% LSS started 24 May 2006, based on probjitter.m
%
ninputspikes = length(oldtimes) ;
newtimes = zeros([1 length(oldtimes)]) ; % set up output array
% generate two sets of random numbers
randarray = rand([1 ninputspikes]) ; % for determining whether a jittered spike will occur
rnarray = rand([1 ninputspikes]) ; % for determining the size of the jitter
outspikeno = 1 ; % index output spikes
for spikeno = 1:ninputspikes
    if randarray(spikeno) < probspike
        % generate an output spike
        deltatime = (rnarray(spikeno) - 0.5) * jitterlength ;        
        newtimes(outspikeno) = oldtimes(spikeno) + deltatime ;
        if (newtimes(outspikeno) > 0) && (newtimes(outspikeno) < duration)
        % disallow if < 0 or beyond end of duration)
            outspikeno = outspikeno + 1 ;
        end
    end
end
if (outspikeno > 1)
    finalnewtimes = newtimes(1:outspikeno - 1) ;
else
    finalnewtimes = [] ;
end

end
        
    
