function f = GPIPcomponent(QC,Phi,w,sigma2,P,finit)
% f: a vector for iteration computation
[A,B] = EffecChanMtxAandB(QC,Phi,sigma2,P);
Abar = Abar_f(A,finit,w,QC);
Bbar = Bbar_f(B,finit,w,QC);
f = PrecoderGPIP(Abar,Bbar,finit);