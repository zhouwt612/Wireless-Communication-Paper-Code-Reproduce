function fty = modulofunc(y,tau)

if nargin < 2
    tau = 2*(abs(3 + 1i*3)+2/2);
end

fty = y - floor((y+tau/2)/tau)*tau;