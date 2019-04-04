function spikesout =  genspikesgaussian(duration, gmean, gstdev, varargin)
% generate a spike train of duration duration, using Gaussian random data.
% we generate normally discributed values 
% with mean gmean and stdev gstdev,both in seconds
%  LSS 8 JUNE 2005
% allow specification of minimum ISI (in seconds)

MinISI = 0.0001 ;
Resolution = 0.00001 ; % 10 microseconds
% Read in arguments. Thanks to David Sterratt for code below!
for i=1:(nargin-3)
    if ischar(varargin{i})
        switch varargin{i}
            case 'MinISI', MinISI = varargin{i+1};
            case 'Resolution', Resolution = varargin{i+1};
        end
    end
end



% estimate length of spikesout
spikesout1 = zeros([1 round(duration * 1.25 / gmean  )]) ;

time = 0 ;
spikeno = 1 ;
while time < duration
    % generate new ISI
    nrand = randn(1) ;
    ISI = (nrand * gstdev) + gmean ;
    if ISI > MinISI % can't have ISI less than some fixed quantity
        spikesout1(spikeno) = time + ISI ;
        time = time + ISI ;
        spikeno = spikeno + 1 ;
    end
end
spikesout = spikesout1(1:spikeno-2 ) ;

end