% BER5=BER4;
% capacity_avg5=capacity_avg4;
% eig_num_avg5=eig_num_avg4;
% round_num_avg5=round_num_avg4;

% load N8K8Kall50Mi1PreSelectSNR9To11.mat BER4 capacity_avg4 eig_num_avg4 round_num_avg4
% 
% BER5=BER4+BER5;
% capacity_avg5=capacity_avg4+capacity_avg5;
% eig_num_avg5=eig_num_avg4+eig_num_avg5;
% round_num_avg5=round_num_avg4+round_num_avg5;

fprintf(1,'%s','ok')

BER4=BER5
capacity_avg4=capacity_avg5;
eig_num_avg4=eig_num_avg5;
round_num_avg4=round_num_avg5;