function ber = invMperM(snr,M,K,datanum,alpha)

no_pow = 1/10^(snr/10);
xdata = zeros(K,datanum);
ydata = zeros(K,datanum);
for ind1 = 1:1:datanum
    x = randsrc(K,1,0:M-1);
    xmod = pskmod(x,M,pi/4);
    h = sqrt(1/2)*(randn(K)+1i*randn(K));
    L = 4; tau = 1;
%     [L,tau] = LtaufinderM(h,xmod,alpha);
    s = chaninv(h,xmod+tau*L);
    noise = sqrt(1/2*no_pow)*(randn(K,1)+1i*randn(K,1));
    snor = s/sqrt(trace(s*s'));
    y = h*snor + noise;
    ymodulo = modulofunc(real(y),tau) + 1i*(modulofunc(imag(y),tau));
    ydemod = pskdemod(ymodulo,M,pi/4);
    xdata(:,ind1) = x;
    ydata(:,ind1) = ydemod;
end
[~,ber] = biterr(xdata,ydata);