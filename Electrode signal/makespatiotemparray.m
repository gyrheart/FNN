function spatiotemparray = makespatiotemparray(L, setorigtimes, setorigvalues, setdifftimes, setdiffvalues, setd2times, setd2values,varargin)
%   produces an array which is 3 * L long, and writes the array to a file
%   if the FileName is supplied
% The parameters are
%   L: the length of the array
%   setorigtimes, setorigvalues: the times and values for the original
%   signal
%   setdifftimes, setdiffvalues: the times and values for the
%   differentialted signal
%   setd2times, setd2values: gthe times and valuyes for the twice
%   differentiated signal
%
%   LSS 23 March 2006.
%
% Note that the output can also be read by PrepareParams

spatiotemparray = zeros([3 L]) ;
% sanity check parameters
if (length(setorigtimes) ~= length(setorigvalues))
    error('original value locations and values not same length') ;
end
% same for differential
if (length(setdifftimes) ~= length(setdiffvalues))
    error('differential value locations and values not same length') ;
end
if (length(setd2times) ~= length(setd2values))
    error('differential value locations and values not same length') ;
end
% now do linear interpolation on all three
spatiotemparray(1,:) = interp1(setorigtimes, setorigvalues, 1:L, 'linear', 0) ;
spatiotemparray(2,:) = interp1(setdifftimes, setdiffvalues, 1:L, 'linear', 0) ;
spatiotemparray(3,:) = interp1(setd2times, setd2values, 1:L, 'linear', 0) ;
%
% do we save the values in a file?
%
for i=1:(nargin-7)
    if ischar(varargin{i})
        switch varargin{i}
            case 'FileName'
                % write the file out
                fhandle = fopen(cell2mat(varargin(i+1)), 'w') ;
                fprintf(fhandle,'%6.4f\t', spatiotemparray) ;
                fclose(fhandle) ;
        end
    end
end