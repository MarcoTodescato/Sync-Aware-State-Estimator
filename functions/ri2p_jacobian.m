function Rx = ri2p_jacobian(x)

[N,M] = size(x);
if N>1 && M> 1
    error('input dimension mismatch')
end

Rx = [];
for i = 1:length(x)
    
    Rx_ii = [cos(angle(x(i))) , -abs(x(i))*sin(angle(x(i))); ...
             sin(angle(x(i))) ,  abs(x(i))*cos(angle(x(i)))];
         
    Rx = blkdiag(Rx,Rx_ii);
    
end

end

