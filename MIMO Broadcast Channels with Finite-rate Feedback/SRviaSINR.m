function sumrate = SRviaSINR(H,G,P)
% sumrate: the total sum rate
% H: channel
% G: precoder
% P: transmit power

[K,M] = size(H);
ratei = zeros(1,K);
for idx1 = 1:1:K
    interference = 0;
    signal_pow = 0;
    for idx2 = 1:1:K
        if idx1 == idx2
            signal_pow = abs(H(idx1,:)*G(:,idx2)/norm(G(:,idx2)))^2;
        else
            interference = interference + abs(H(idx1,:)*G(:,idx2)/norm(G(:,idx2)))^2;
        end
    end
    SINRi = (P/M)*signal_pow/(1+(P/M)*interference);
    ratei(1,idx1) = log2(1+SINRi);
end
sumrate = sum(ratei);