function sumrate = one_transmission(K,M,snr,B,Hnum)
% SR: sum rate for one channel realization
% K: the number of UE
% M: the number of transmit antennas
% snr: signal-to-noise ratio
% B: feedback bits

% frame_len: the length of one frame
% In this simulation, assuming the power of noise is 1.
sigma2 = 1;
snr_linear = 10^(snr/10);
P = sigma2*snr_linear;
C = RVQcodebook(M,B);
sumratedata = zeros(1,Hnum);
for idx1 = 1:1:Hnum
%     C = RVQcodebook(M,B);
    H = channel(K,M);
    QuantH = quantizedchannel(H,C);
    Pre = RegZFprecoder(QuantH,P);
%     Pre = QuantH'/(QuantH*QuantH');
%     sumratedata(1,idx1) = SRviaSINRnorm(H,Pre,P);        
    sumratedata(1,idx1) = SRergodic(H,Pre,P);
end
% if sum(isnan(sumratedata(1,:))) == 0
%     
sumrate = M*mean(sumratedata);