function f = PrecoderGPIP(Abar,Bbar,finit)
% f: (NK, 1) computed precoder through GPIP algorithm
% Abar: (NK, NK) a matrix for computing f
% Bbar: (NK, NK) a matrix for computing f
% finit: (NK, 1) initial f
fb = inv(Bbar)*Abar*finit;
f = fb/norm(fb,2);