function u = EMMSE(MT,R,K,KK,H,sigma2,FB_bits)

err_UpperBound=2^(-FB_bits/(MT-1));

% 根据调度结果，读入被调度者的信道系数
H_now=zeros(K,MT);
for i1=1:K,
    H_now(i1,:)=H((KK(1,i1)-1)*R+1,:);
end

% pre-coding
% The transmitter only knows the quantized channel correlation matrix.
u=zeros(MT,K);      % single stream
u=inv(H_now'*H_now+((1+(K/sigma2)*err_UpperBound)/((K/sigma2)*(1-err_UpperBound)))*eye(MT))*H_now';
u=u/(sqrt(trace(u'*u)));
