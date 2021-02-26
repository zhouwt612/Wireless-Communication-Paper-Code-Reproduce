clc
clear
close

t = 3;
r = 3;
N = 3;
Pow = 1;
snr = 5;

P = orth_pilot(N,t,Pow);
H = sqrt(1/2)*(randn(r,t)+1i*randn(r,t))
V = sqrt(1/2*(Pow/10^(snr/10)))*(randn(r,t)+1i*randn(r,t));
S = H*P+V;
Hls = LSestimator(S,P)