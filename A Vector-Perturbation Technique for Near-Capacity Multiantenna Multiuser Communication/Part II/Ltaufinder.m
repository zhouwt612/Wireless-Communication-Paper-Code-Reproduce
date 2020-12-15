function Lopt = Ltaufinder(H,u,K,Lnum,tau)
% To find L through Eq. 8

if nargin < 5
    tau = 2*(abs(3 + 1i*3)+2/2);
end

L = round(rand(K,Lnum)*20)+1i*round(rand(K,Lnum)*20);
uHHudata = zeros(K,Lnum);
for idx = 1:Lnum
    uHHu = (u+tau*L(:,idx))'/(H*H')*(u+tau*L(:,idx));
    uHHudata(:,idx) = uHHu;
end
[~,I] = min(uHHudata(1,:));
Lopt = L(:,I);