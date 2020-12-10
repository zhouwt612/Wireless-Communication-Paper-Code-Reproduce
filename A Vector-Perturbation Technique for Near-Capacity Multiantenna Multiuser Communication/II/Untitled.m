clc
clear
close

K = 10;
H = sqrt(1/2)*(randn(K)+1i*randn(K));
invHH = inv(H*H');

tau = 2*(abs(3+1i*3)+1);
l_prime = 1:1:100;
x = randsrc(K,1,0:15);
xmod = qammod(x,16);

