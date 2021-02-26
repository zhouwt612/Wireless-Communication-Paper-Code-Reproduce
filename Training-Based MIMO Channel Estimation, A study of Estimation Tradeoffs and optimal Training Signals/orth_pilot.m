function P = orth_pilot(N,t,pow)
% Generate an orthogonal pilot matrix
% N: the length of pilots
% t: the number of transmit antennas
% pow: power constraint

P = zeros(t,N);
Wn = exp(1i*2*pi/N);

for idx1 = 1:1:t
    for idx2 =1:1:N
        P(idx1,idx2) = Wn^((idx1-1)*(idx2-1));
    end
end
P = sqrt(pow/N/t)*P;