clc
clear
p0 = 1;

for i = 1:1000
    disp([num2str(i) ' realizations out of ' num2str(100)]);
    h = sqrt(1/2)*(randn(2,2)+1i*randn(2,2));
    R = eye(2);
    [v,lambda] = eig(h'/R*h);
    [d,ind] = sort(diag(lambda),'descend');
    LAMBDA = lambda(ind,ind);
    V = v(:,ind);
    V2x2 = V(:,[1 2]);
    LAMBDA2x2 = LAMBDA([1 2],[1 2]);
    % precoder
    w = eye(2);
    mu12 = trace(LAMBDA2x2^(-1/2)*w^(1/2))/(p0+trace(LAMBDA2x2^(-1)));
    phi_f = mu12^(-1)*LAMBDA2x2^(-1/2)*w^(1/2) - LAMBDA2x2^(-1)
    %     if phi_f(4) < 0
    %         phi_f((phi_f<0)) = 0;
    %         phi_f = phi_f/phi_f(1);
    %     end
    phi_f((phi_f<0)) = 0;
    phi_f = phi_f^(1/2);
    %     rho = LAMBDA2x2*w
    %     mu = mu12^2
    
    
    
    F = V2x2 * phi_f;
    
    Tphif = trace(phi_f*phi_f')
end

