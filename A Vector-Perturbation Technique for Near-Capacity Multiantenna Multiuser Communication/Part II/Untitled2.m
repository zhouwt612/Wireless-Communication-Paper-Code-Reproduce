clc
clear
% close

M = 10;
K = 10;

H = sqrt(1/2)*(randn(M,K)+1i*randn(M,K));

tau = 2*(abs(3 + 1i*3)+2/2);
x = randsrc(K,1,0:16-1);
u = qammod(x,16);

L_len = 20000;
L = round(rand(K,L_len)*50)+1i*round(rand(K,L_len)*50);

uHHudata = zeros(K,L_len);
for idx = 1:L_len
    uHHu = (u+tau*L(:,idx))'/(H*H')*(u+tau*L(:,idx));
    uHHudata(:,idx) = uHHu;
end

[M,I] = min(uHHudata(1,:))
L(:,I)