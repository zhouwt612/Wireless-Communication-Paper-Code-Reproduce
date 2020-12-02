function G = RxWODE(H,P,R_s,R_n)

G = R_s'*P'*H'*(H*P*R_s'*P'*H'+R_n')^(-1);