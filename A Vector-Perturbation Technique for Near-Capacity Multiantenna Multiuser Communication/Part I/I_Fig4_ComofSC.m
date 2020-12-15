clc
clear
close

syms t
y = @(t)(exp(-t)./t);

K = 2:2:16;
rho = 10;

hnum = 10000;
Csummean = zeros(1,length(K));
Cregdata = zeros(1,length(K));
Ccidata = zeros(1,length(K));
for ind1 = 1:1:length(K)
    disp(['The number of users: ' num2str(K(ind1))])
    Csumspa = zeros(1,hnum);
    Cregspa = zeros(1,hnum);
    for ind2 = 1:1:hnum
        % sum capacity
        h = sqrt(1/2)*(randn(K(ind1))+1i*randn(K(ind1)));
        Ik = eye(K(ind1));
        D = (1/K(ind1))*Ik;
        Csum = log2(det(Ik+rho*h'*D*h));
        Csumspa(1,ind2) = Csum(:);
        % regularized inversion
        SINR = SINRcal(h,K(ind1),K(ind1)/rho,rho);
        Creg = K(ind1)*log2(1+SINR);
        Cregspa(1,ind2) = Creg;
        % channel inversion
        E1Krho = integral(y,K(ind1)/rho,Inf);        
        Cci = K(ind1)*exp(K(ind1)/rho)*E1Krho;
    end
    Csummean(1,ind1) = mean(Csumspa);
    Cregdata(1,ind1) = mean(Cregspa);
    Ccidata(1,ind1) = Cci;
end

Csummean
Cregdata
Ccidata

plot(K,real(Csummean),'r--x','linewidth',1.5)
hold on
plot(K,Cregdata,'g','linewidth',1.5)
plot(K,Ccidata,'b-.s','linewidth',1.5)
xlabel('K');ylabel('Capacity')
legend('Sum Capacity','Regularized Inversion','Channel Inversion')
% axis([2 16 0 45])