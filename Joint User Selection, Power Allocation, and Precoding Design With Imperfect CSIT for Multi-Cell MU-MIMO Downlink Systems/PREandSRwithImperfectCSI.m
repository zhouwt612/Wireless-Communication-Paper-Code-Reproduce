function [Bwmmse,SR] = PREandSRwithImperfectCSI(H,F,Etx,u,P,Q,K,itenum)
B = F'*sqrt(Etx/trace(F'*F));
Rn = nointercov(H,B,P,Q,K);
% WSRBF-WMMSE
Ammse = zeros(Q*K,Q*K);
W = zeros(Q*K,Q*K);
Datarate = zeros(1,itenum);
for index = 1:1:itenum
    
    % Update Ammse
    for i = 1:1:K
        Ammse(((i-1)*Q+1):i*Q,((i-1)*Q+1):i*Q) = RxMMSE(F(((i-1)*Q+1):i*Q,:),...
            B(:,((i-1)*Q+1):i*Q),Rn(:,:,i));
    end
    % Update W
    for j = 1:1:K
        W(((j-1)*Q+1):j*Q,((j-1)*Q+1):j*Q) = MSEweight(F(((j-1)*Q+1):j*Q,:)...
            ,B(:,((j-1)*Q+1):j*Q),Rn(:,:,j));
    end
    % Update Bmmse
    B = TxWMMSE(F,Ammse,W,Etx);
    
    % Compute the sum rate
    Rn = nointercov(H,B,P,Q,K);
    sumrate = SUMrateMISO(H,B,Rn,u,K);
    Datarate(index) = sumrate;
end
SR = Datarate;
Bwmmse = B;