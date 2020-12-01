function total_vector_error = computeTVE(vecA, vecB)
% Total vector error between two vectors
% 
% Args:
%   vecA: first vector
%   vecB: second vector
%
% Returns:
%   TVE

total_vector_error = mean(mean(abs(vecA - vecB) ./ abs(vecB)));

end

