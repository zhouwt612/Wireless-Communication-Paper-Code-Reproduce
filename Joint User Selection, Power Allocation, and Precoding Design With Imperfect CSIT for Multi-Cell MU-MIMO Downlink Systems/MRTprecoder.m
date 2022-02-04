function precoder = MRTprecoder(H,P)
PRE = H';
b = sqrt(P/trace(PRE*PRE'));
precoder = b*PRE;