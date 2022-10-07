clc
clear
close

P = 4;
Q = 2;
K = 2;

SNR = 10;
H = sqrt(10^(SNR/10)/2)*(randn(Q*K,P)+1i*randn(Q*K,P));
B = randn(P,Q*K)+1i*randn(P,Q*K);
% B = ones(P,Q*K);
u = ones(1,K);

BD = B;
Etx = 1;

% WSRBF-WMMSE
Ammse = zeros(Q*K,Q*K);
W = zeros(Q*K,Q*K);
Datarate = zeros(1,30);
iteration = 1:1:30;
for index = iteration
    Rn = nointercov(H,B,P,Q,K);
    % Update Ammse
    for i = 1:1:K
        Ammse(((i-1)*Q+1):i*Q,((i-1)*Q+1):i*Q) = RxMMSE(H(((i-1)*Q+1):i*Q,:),...
            B(:,((i-1)*Q+1):i*Q),Rn(:,:,i));
    end
    % Update W
    for j = 1:1:K
        W(((j-1)*Q+1):j*Q,((j-1)*Q+1):j*Q) = MSEweight(H(((j-1)*Q+1):j*Q,:)...
            ,B(:,((j-1)*Q+1):j*Q),Rn(:,:,j),u(j));
    end
    % Update Bmmse
    B = TxWMMSE(H,Ammse,W,Etx);
    
    % Compute the sum rate
    sumrate = SUMrate(H,B,Rn,u,Q,K);
    Datarate(index) = sumrate;
end

% WSRBF-WMMSE-D
DatarateD = zeros(1,30);
for indexD = iteration
    RnD = nointercov(H,BD,P,Q,K);
    % Update Ammse
    for iD = 1:1:K
        AmmseD(((iD-1)*Q+1):iD*Q,((iD-1)*Q+1):iD*Q) = RxMMSE(H(((iD-1)*Q+1):iD*Q,:),...
            BD(:,((iD-1)*Q+1):iD*Q),RnD(:,:,iD));
    end
    % Update W
    for jD = 1:1:K
        WD(((jD-1)*Q+1):jD*Q,((jD-1)*Q+1):jD*Q) = MSEweightD(H(((jD-1)*Q+1):jD*Q,:)...
            ,BD(:,((jD-1)*Q+1):jD*Q),RnD(:,:,jD),u(j));
    end
    % Update Bmmse
    BD = TxWMMSE(H,AmmseD,WD,Etx);
    
    % Compute the sum rate
    sumrateD = SUMrate(H,BD,RnD,u,Q,K);
    DatarateD(indexD) = sumrateD;
end

figure
plot(iteration,Datarate,'r-*',iteration,DatarateD,'g-o')
xlabel('Iteration number n');ylabel('Sum-rate');
legend('WSRBF-WMMSE','WSRBF-WMMSE-D');
grid on
