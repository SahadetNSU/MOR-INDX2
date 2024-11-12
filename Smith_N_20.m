%% projection free solution initialization 


tic 

iter1 = 8;      % total iteration number (30/50 for 20*20, 100/200 for 1200)

% c_mat = [E1 A2;A3 A4];
Zt_r = [E1 A2;A3 A4]\[B1; B2];
Z_r = Zt_r(1:size(E1,1), 1:end);
sigma_r = Zt_r(size(E1,1)+1:end,1:end);

Zt_l = [E1' A2;A3 A4]\[C1'; C2'];
Z_l = Zt_l(1:size(E1,1), 1:end);
sigma_l = Zt_l(size(E1,1)+1:end,1:end);

R = cell(1,iter1);
R{1} = Z_r;

L = cell(1,iter1);
L{1} = Z_l;

disp('normalized residual norms for controllability and observability Gramians')
res1_cc(1)=norm(E_m*R{1}*R{1}'*E_m' - A_m*R{1}*R{1}'*A_m' - B_m*B_m','fro')/...
    norm(B_m*B_m','fro');
fprintf(1,'step: 1 normalized residual: %d\n',res1_cc(1));

res1_oc(1)=norm(E_m'*L{1}*L{1}'*E_m-A_m'*L{1}*L{1}'*A_m- C_m'*C_m,'fro')/...
    norm(C_m'*C_m,'fro');
fprintf(1,'step: 1  normalized residual: %d\n',res1_oc(1));


for i = 2:iter1
    
    Z_new_r = A1*Z_r;
    Zt_r = [E1 A2;A3 A4]\[Z_new_r; B2];
    Z_r = Zt_r(1:size(A1,1), 1:end);
    
    Z_new_l = A1'*Z_l;
    Zt_l = [E1' A2;A3 A4]\[Z_new_l; C2'];
    Z_l = Zt_l(1:size(E1,1), 1:end);
    
    R{i} = [R{i-1} Z_r];
    L{i} = [L{i-1} Z_l];
    
    res1_cc(i)=norm(E_m*R{i}*R{i}'*E_m' - A_m*R{i}*R{i}'*A_m' - B_m*B_m','fro')/...
    norm(B_m*B_m','fro');
    fprintf(1,'step: %4d  normalized residual: %d\n',i,res1_cc(i));
    
    res1_oc(i)=norm(E_m'*L{i}*L{i}'*E_m-A_m'*L{i}*L{i}'*A_m- C_m'*C_m,'fro')/...
    norm(C_m'*C_m,'fro');
    fprintf(1,'step: %4d  normalized residual: %d\n',i,res1_oc(i));
end

toc

clear Zt Z sigma Z_new ;
