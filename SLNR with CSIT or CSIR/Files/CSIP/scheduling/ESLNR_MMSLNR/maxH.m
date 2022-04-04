function KK = maxH(N,Mi,S,Kall,K,H)

channel_norm=zeros(1,Kall);
HHH=zeros(Mi,N);
for i1=1:Kall,
    HHH=H((i1-1)*Mi+1:i1*Mi,:);
    channel_norm(1,i1)=trace(HHH'*HHH);
end

KK=zeros(1,K);

[ch_norm,index]=another_sort(channel_norm);

for i1=1:K,
    KK(1,i1)=index(1,i1);
end