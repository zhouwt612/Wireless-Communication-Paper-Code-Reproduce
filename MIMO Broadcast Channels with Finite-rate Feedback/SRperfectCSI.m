function sumrate = SRperfectCSI(M,P,c)
% M: the number of transmit antennas
% P: transmit power constraint
% c: a constant depending on the channel realization H

if nargin < 3
    c = 0;
end

sumrate = M*log(P) + c;