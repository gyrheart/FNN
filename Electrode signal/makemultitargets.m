function multitargets =  makemultitargets(original, n_sets)
% produce a number of targets by linearly interpolating between the
% original pair
% n_sets > 1 (pointless otherwise)
[t1 t2 t3] = size(original) ;
% n.b. t1 should be 2: if not the 1st 2 are used
% t2 must be 3
multitargets = zeros([n_sets 2 t2 t3]) ;
delta = original(2,:,:) - original(1,:,:) ;
increment = delta/n_sets ;
for i=1:n_sets
    
    multitargets(i,2,:,:) = original(2,:,:) ;
    multitargets(i,1,:,:) = original(1,:,:) + (i - 1) * increment ;
end