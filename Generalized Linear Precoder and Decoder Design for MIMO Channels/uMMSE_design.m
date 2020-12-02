function [F,G] = uMMSE_design(h,p0,B,sigmasq,w)
% compute the precoder and encoder of the channel through unweighted MMSE

R = sigmasq*eye(B);
w = 1/w*eye(B);
[v,lambda] = eig(h'/R*h);
[d,ind] = sort(diag(lambda),'descend');
LAMBDA = diag(d);
V = v(:,ind);
V2x2 = V(:,1:B);
LAMBDA2x2 = LAMBDA(1:B,1:B);

% precoder
mu12 = trace(LAMBDA2x2^(-1/2)*w^(1/2))/(p0+trace(LAMBDA2x2^(-1)));
phi_f = mu12^(-1)*LAMBDA2x2^(-1/2)*w^(1/2) - LAMBDA2x2^(-1);
phi_f((phi_f<0)) = 0;
phi_f = phi_f^(1/2);
F = V2x2 * phi_f;

%decoder
phi_g = mu12*LAMBDA2x2^(-1/2)*w^(-1/2)-(mu12^2)*LAMBDA2x2^(-1)*w^(-1);
phi_g((phi_g<0)) = 0;
phi_g = phi_g^(1/2)*LAMBDA2x2^(-1/2);
G = phi_g*V2x2'*h'/R;
