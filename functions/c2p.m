function Z = c2p(z)
% COMPLE2REIM given a complex valued number, vector or matrix, gives back
% its polar representation as matrix.
%
% Example:
%
% z = a + 1i*b ---> Z = [abs(z); ang(z)];

% dimensions
[NumberOfRow , NumberOfColumn, NumberOfSlice] = size(z);

% init the final matrix
Z = zeros(2*NumberOfRow , NumberOfColumn, NumberOfSlice);

% loop over the rows and columns to fill the element of Z
for s = 1:NumberOfSlice
    for r = 1:NumberOfRow
        for c = 1:NumberOfColumn
            
            % real part of the selected element
            abs_z_ij = abs(z(r,c,s));
            
            % imaginary part of the selected element
            ang_z_ij = angle(z(r,c,s));
            
            % fill the 2x2 sub-block of Z with the real and imag part of z_ij
            Z(2*r-1:2*r,c,s) = [abs_z_ij ; ang_z_ij];
            
        end
    end
end
