function ber = reginvtransmission(snr,K,datanum)

no_pow = 1/10^(snr/10);
xdata = zeros(K,datanum);
ydata = zeros(K,datanum);
for ind1 = 1:1:datanum
    x = randsrc(K,1,0:3);
    xmod = qammod(x,4);
    h = sqrt(1/2)*(randn(K)+1i*randn(K));
    s = regchaninv(h,xmod,K*no_pow,K);
    noise = sqrt(1/2*no_pow)*(randn(K,1)+1i*randn(K,1));
    snor = s/sqrt(trace(s*s'));
    y = h*snor + noise;
    ydemod = qamdemod(y,4);
    xdata(:,ind1) = x;
    ydata(:,ind1) = ydemod;
end
[~,ber] = biterr(xdata,ydata);