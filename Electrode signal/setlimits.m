function sigout = setlimits(samplevalues, MinVal, MaxVal)
% linearly scale the values in samplevalues to b between MinVal and MaxVal
% LSS 13 June 2005
%
origmin = min(samplevalues) ;
origmax = max(samplevalues) ;
% linear stretch: multiply by k = (MaxVal-MinVal)/(origmax - origmin) then
% add MinVal - k * origmin
if ((origmax - origmin) == 0)
    sigout = zeros(size(samplevalues)) ;
else
    stretchk = (MaxVal-MinVal)/(origmax - origmin) ;
    sigout1 = samplevalues * stretchk ;
    sigout = sigout1 + (MinVal - stretchk * origmin) ;
end