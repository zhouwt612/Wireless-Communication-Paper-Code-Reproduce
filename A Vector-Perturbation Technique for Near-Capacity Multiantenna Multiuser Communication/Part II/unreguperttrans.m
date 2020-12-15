function ber = unreguperttrans(K,M,Modu,snr,datanum,Lnum,tau)

if nargin < 6
    tau = 2*(abs(3 + 1i*3)+2/2);
end

no_pow = 1/10^(snr/10);
xdata = zeros(K,datanum);
ydata = zeros(K,datanum);
for idx = 1:1:datanum
    x = randsrc(K,1,0:Modu-1);
    u = qammod(x,Modu);
    H = sqrt(1/2)*(randn(M,K)+1i*randn(M,K));
    Lopt = Ltaufinder(H,u,K,Lnum,tau);
    s = chaninv(H,u+tau*Lopt);
    noise = sqrt(1/2*no_pow)*(randn(K,1)+1i*randn(K,1));
    snor = s/sqrt(trace(s*s'));
    y = H*snor + noise;
    yI = modulofunc(real(y),tau);
    yQ = modulofunc(imag(y),tau);
    ymodulo = yI + 1i*yQ;
    ydemod = qamdemod(ymodulo,Modu);
    xdata(:,idx) = x;
    ydata(:,idx) = ydemod;
end
[~,ber] = biterr(xdata,ydata);