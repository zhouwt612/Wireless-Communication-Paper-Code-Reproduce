function ber = reginvM(snr,M,K,datanum)

no_pow = 1/10^(snr/10);
xdata = zeros(K,datanum);
ydata = zeros(K,datanum);
for ind1 = 1:1:datanum
    x = randsrc(K,1,0:M-1);
    xmod = qammod(x,M);
    h = sqrt(1/2)*(randn(K)+1i*randn(K));
    s = regchaninv(h,xmod,K*no_pow,K);
    noise = sqrt(1/2*no_pow)*(randn(K,1)+1i*randn(K,1));
    snor = s/sqrt(trace(s*s'));
    y = h*snor + noise;
    ydemod = qamdemod(y,M);
    xdata(:,ind1) = x;
    ydata(:,ind1) = ydemod;
end
[~,ber] = biterr(xdata,ydata);