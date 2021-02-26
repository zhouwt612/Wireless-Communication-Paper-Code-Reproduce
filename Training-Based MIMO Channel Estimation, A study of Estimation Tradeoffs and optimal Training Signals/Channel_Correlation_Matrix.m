clc
% clear
% close

hnum = 10000;
r = 3;
t = 3;

h = sqrt(1/2)*(randn(r,t,hnum)+1i*randn(r,t,hnum));

HHH = zeros(t,t,hnum);
for idx1 = 1:1:hnum
    HHH(:,:,idx1) = h(:,:,idx1)'*h(:,:,idx1);
end

Rh = zeros(t,t);
for idx2 = 1:1:t
    for idx3 = 1:1:t
        Rh(idx2,idx3) = mean(HHH(idx2,idx3,:));
    end
end
Rh