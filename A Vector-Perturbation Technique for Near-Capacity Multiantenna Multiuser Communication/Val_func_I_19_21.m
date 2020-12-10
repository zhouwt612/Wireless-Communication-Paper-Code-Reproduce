clc
clear
close

% function (19) (20) (21)
H = sqrt(1/2)*(randn(4)+1i*randn(4));
HH = H*H';
[Q1,lambda] = eig(HH);
HH_alpha = H*H'/(H*H'+eye(4))
HH_alpha2 = Q1*(lambda/(lambda+eye(4)))*Q1'

HH_alpha - HH_alpha2
