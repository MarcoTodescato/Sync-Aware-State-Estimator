function Z = c2ri(z, varargin)
% COMPLE2REIM given a complex valued number, vector or matrix, gives back
% its real and imag representation as matrix.
%
% Example:
%
% z = a + 1i*b ---> Z = [a -b ; b  a];

% dimensions
[NumberOfRow , NumberOfColumn, NumberOfSlice] = size(z);

% init the final matrix
Z = zeros(2*NumberOfRow , 2*NumberOfColumn, NumberOfSlice);

% loop over the rows and columns to fill the element of Z
for s = 1:NumberOfSlice
    for r = 1:NumberOfRow
        for c = 1:NumberOfColumn
            
            % real part of the selected element
            re_z_ij = real(z(r,c,s));
            
            % imaginary part of the selected element
            im_z_ij = imag(z(r,c,s));
            
            % fill the 2x2 sub-block of Z with the real and imag part of z_ij
            Z(2*r-1:2*r,2*c-1:2*c,s) = [re_z_ij , -im_z_ij ; im_z_ij , re_z_ij];
            
        end
    end
end

if nargin > 1
    flag = varargin{:};
    
    switch flag
        case 'v'
            Z = Z(:,1,:);
            
        case 'evol'
            Z = Z(:,1:2:end,:);
            
        otherwise
            
    end
end

