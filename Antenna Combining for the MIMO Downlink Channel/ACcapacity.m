function C = ACcapacity(H,Pre,pow)
% The capacity with equal power loading
% H: the effective channel (size: M,K)
% Pre: precoder, (size: M,K)
% pow: transmit power

[M,K] = size(H);
R = zeros(1,K);
for idx1 = 1:1:K
    powi = 0;
    powineqk = 0;
    for idx2 = 1:1:K
        if idx1 == idx2
            powi = (pow/M)*abs(H(:,idx1)'*Pre(:,idx1))^2;
        else
            powineqk = powineqk + (pow/M)*abs(H(:,idx1)'*Pre(:,idx2))^2;
        end
    end
    R(idx1) = log2(1 + powi/powineqk);
end
C = sum(R);