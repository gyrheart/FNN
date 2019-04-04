function correlateddata =  correlated_temporal_mac(n_correlateds, timepoints, dirstring)
% function which initialises the time envelopes for the correlated neurons.
% There are n_correlateds correlated neurons, and timepoints points in time. The
% actual timing is set in the calling function.
% LSS started 10 3 2006
%
correlateddata = zeros([n_correlateds 3 timepoints]) ;

% LSS started 10 3 2006
%

%

% open the file
for i = 1: n_correlateds
    fno = i - 1 ;
    fname = [dirstring '/correlated_temporal_' int2str(fno) '.dat'] ;
    fhandle = fopen(fname, 'r') ; % open for reading
    if fhandle == -1
        disp(['File ' dirstring fname ' cannot be opened']) ;
        return ;
    end
    correlateddata(i,:,:) = fscanf(fhandle,'%f', [3 timepoints]) ; 
    fclose(fhandle) ;
end