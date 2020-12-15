clc
clear
close

K = 4:2:20;
cir_num = 5000;
Elogdata = zeros(4,length(K));
for ind1 = 1:1:length(K)
    Edata = zeros(4,cir_num);
    for ind2 = 1:1:cir_num
        h = sqrt(1/2)*(randn(K(ind1))+1i*randn(K(ind1)));
        hhH = inv(h*h');
        [~,E] = eig(hhH);
        Ediag = diag(E);
        Ediag = Ediag(end-3:end,:);
        Edata(:,ind2) = Ediag;
    end
    Em = mean(Edata,2);
    Elog = log10((1/K(ind1))*Em);
    Elogdata(:,ind1) = Elog;
end
Elogdata

figure
plot(K,Elogdata(1,:),'linewidth',2);hold on;
plot(K,Elogdata(2,:),'linewidth',2);plot(K,Elogdata(3,:),'linewidth',2)
plot(K,Elogdata(4,:),'linewidth',2)
xlabel('Dimension (K)');ylabel('Log10[(1/K)Eig]');axis([4 20 -2 2])
title('Four largest Eigenvalues of inv(HH^T)')