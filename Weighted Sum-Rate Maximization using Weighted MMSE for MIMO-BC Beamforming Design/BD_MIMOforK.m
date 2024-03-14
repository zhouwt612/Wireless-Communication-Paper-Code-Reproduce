function P = BD_MIMOforK(H,pow)
% H: channel realization (size: (M,N,K))
% precoder: precoder matrix (size: (M,N,K))

[M,N,K] = size(H);
P = zeros(M,N,K);
for idx1 = 1:1:K
    Hstack = [];
    for idx2 = 1:1:K
        if idx1 ~= idx2
            Hstack = [Hstack H(:,:,idx2)];
        end
    end
    precoder = null(Hstack');
%     [U S V] = svd(Hstack');
%     precoder = V(:,end-N+1:end);
    P(:,:,idx1) = sqrt(pow/K/trace(precoder*precoder'))*precoder;
end