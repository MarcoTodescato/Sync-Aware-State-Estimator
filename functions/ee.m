function v = ee(indices, l)
% v = ee(indices, n)
%  returns a n-dimensional vector of 1 and -1 in the positions defined in 
%  indices. Negative indices means negative entries.
%  
%  EXAMPLE
%  ee([2 -3 -5], 7) = [0 1 -1 0 -1 0 0]'
% 

  v = zeros(l,1);
  v(indices(indices > 0)) = 1;
  v(-indices(indices < 0)) = v(-indices(indices < 0)) -1;

end
