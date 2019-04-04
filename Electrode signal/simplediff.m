function diffout = simplediff(sigin, smoothlength)
%
% performs a simple differential operation. 
% smooths input  sigin first
% uses  hamming window (should be selectable) of length 2* smoothlength
% return vector is same length as input vector: always 0 for first and last
% smoothlength elements
%
% LSS 10 6 2005.
%
% attempting optimisation (8 2 2006)

triangwindow = hamming(smoothlength * 2) ;
convsigin = conv(sigin, triangwindow) ;
% would initialising diffout help? Yes. 
diffout = zeros([1 length(sigin)]) ;
for sampleno = 1:smoothlength
    diffout(sampleno) = 0 ;
end
for sampleno = 2 * smoothlength + 1:length(convsigin)-  2 * smoothlength 
    diffout(sampleno - smoothlength) = convsigin(sampleno+1) - convsigin(sampleno) ;
end
for sampleno = length(convsigin)-  2 * smoothlength + 1:length(convsigin)-   smoothlength + 1
    diffout(sampleno - smoothlength) = 0 ;
end