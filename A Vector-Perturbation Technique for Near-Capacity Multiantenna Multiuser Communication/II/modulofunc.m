function fty = modulofunc(y,tau)
% y: data  tau: a positive real value

fty = y - floor((y+tau/2)/tau)*tau;