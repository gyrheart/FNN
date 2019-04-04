function finalnewtimes = probjitter(oldtimes, probspike, jitterstd, duration)
%
% produces a new list of times from the old list. Each output spike is
% produced from each input spike with probability probspike, at a time
% determined from a gaussian distribution with mean 0 and std jitterstd
%
% no additional spikes cuirrently produced.
%
% LSS started 10 June 2005
%
ninputspikes = length(oldtimes) ;
newtimes = zeros([1 length(oldtimes)]) ; % set up output array
% generate a set of random numbers
randarray = rand([1 ninputspikes]) ;
rnarray = randn([1 ninputspikes]) ;
outspikeno = 1 ; % index output spikes
for spikeno = 1:ninputspikes
    if randarray(spikeno) < probspike
        % generate an output spike
        deltatime = rnarray(spikeno) * jitterstd ;        
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
        
    
