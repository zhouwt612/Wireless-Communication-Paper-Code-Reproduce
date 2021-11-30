function sumrate1 = SRergodic(H,Pre)
% sumrate: the total sum rate
% H: channel
% G: precoder
% P: transmit power

[K,M] = size(H);
interference = 0;
i = 1;
for idx = 1:1:K
    if idx == i
        sigpow = abs(H(i,:)*Pre(:,idx))^2;
    else
        interference = interference + abs(H(i,:)*Pre(:,idx))^2;
    end
end
sumrate1 = log2(1+sigpow/(1+interference));