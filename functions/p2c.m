function z = p2c(Z)
% COMPLE2REIM given a vector or matrix in polar coordinate gives back
% its complex representation.
%
% Example:
%
% Z = [a ; b ] ---> z = a*exp(1i*b);

% dimensions
[NumberOfRow , NumberOfColumn, NumberOfSlice] = size(Z);

% init the final matrix
z = zeros(NumberOfRow/2 , NumberOfColumn, NumberOfSlice);

% loop over the rows and columns to fill the element of Z
for s = 1:NumberOfSlice
    for r = 1:NumberOfRow/2
        for c = 1:NumberOfColumn
            
            z(r,c,s) = Z(2*r-1,c,s)*exp(1i*Z(2*r,c,s));
            
        end
    end
end
