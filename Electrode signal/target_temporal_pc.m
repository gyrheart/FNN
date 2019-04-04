function targetdata =  target_temporal_mac(n_targets, timepoints, dirstring)
% function which initialises the time envelopes for the target neurons.
% There are n_targets target neurons, and timepoints poionts in time. The
% actual timing is set in the calling function.
%
% This version simply reads a file called target_temporal_i.dat, (where i 
% goes from 0 to n_targets - 1) created
% using the PrepareParams Objctive C program (Macintosh only)
%
% LSS started 10 3 2006
%
targetdata = zeros([n_targets 3 timepoints]) ;
%

% open the file
for i = 1: n_targets
    fno = i - 1 ;
    fname = [dirstring '/target_temporal_' int2str(fno) '.dat'] ;
    fhandle = fopen(fname, 'r') ; % open for reading
    if fhandle == -1
        disp(['File ' dirstring fname ' cannot be opened']) ;
        return ;
    end
    targetdata(i,:,:) = fscanf(fhandle,'%f', [3 timepoints]) ; 
end
  
