function [H,Hhat,E] = channel(M,N,K,sigmae)
% M: the number of transmit antennas
% N: the number of receive antennas
% K: the number of users
% sigmae: variance of estimation errors

Hhat = sqrt(1-sigmae)*(randn(M,N,K) + 1i*randn(M,N,K))/sqrt(2);
E = sqrt(sigmae)*(randn(M,N,K) + 1i*randn(M,N,K))/sqrt(2);
H = Hhat + E;