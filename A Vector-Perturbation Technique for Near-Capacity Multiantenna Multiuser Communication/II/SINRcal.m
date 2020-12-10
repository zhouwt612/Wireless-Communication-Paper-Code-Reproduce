function SINR = SINRcal(h,K,alpha,rho)

hhH = h*h';
[~,S] = eig(hhH);
sigma2 = 1/rho;

%E_gamma
E_gamma = trace(S/((S+alpha*eye(K))^2));
De_top = trace(S/(S+alpha*eye(K)))^2;
Unde_bottom = K*trace((S/(S+alpha*eye(K)))^2)-trace(S/(S+alpha*eye(K)))^2;
No_bottom = sigma2*(K^2)*E_gamma;

SINR = De_top/(No_bottom+Unde_bottom);


