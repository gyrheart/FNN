function spikesout =  genspikespoisson(duration, lambda, meanISI, varargin)
% generate a spike train of duration duration, using Poison random data
% Because the dispersion depends on lambda, we'll set lambda and the mean
% ISI independently. 
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

multiplier = meanISI/lambda ; % use to adjust ISI times

% estimate length of spikesout
spikesout1 = zeros([1 round(duration * 1.25 / meanISI  )]) ;

time = 0 ;
spikeno = 1 ;
while time < duration
    % generate new ISI
    prand = poissrnd(lambda) ;
    ISI = prand * multiplier ;
    if ISI > MinISI % can't have ISI less than some fixed quantity
        spikesout1(spikeno) = time + ISI ;
        time = time + ISI ;
        spikeno = spikeno + 1 ;
    end
end
spikesout = spikesout1(1:spikeno-2 ) ;

end
    