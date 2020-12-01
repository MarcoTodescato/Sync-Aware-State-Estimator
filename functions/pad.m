function v = pad(n, pos, val)
% pad(n,pos,val)
%   returns a n-dimensional array with the elements of val in position pos,
%   and zeros elsewhere

v = zeros(n,1);
v(pos) = val;
