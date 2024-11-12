
%sys_proj =  dss(A_n,B_n,C_n,0,E_n,1);
sys_ori = dss(A,B,C,0,E,1);
sys_r = dss(Ar1,Br1,Cr1,0,Er1,1);


%A_cyc = A_n; B_cyc = B_n; C_cyc = C_n; E_cyc = E_n;
A_cyc_ori = A; B_cyc_ori = B; C_cyc_ori = C; E_cyc_ori = E;
A_r = Ar1; B_r = Br1; C_r = Cr1; E_r = Er1;


points = 501;
freqlow = 10^(2);
freqhig = 10^(8);
%freqhig =100000
%linear space 
linst = (freqhig-freqlow)/(points-1);
freq = freqlow:linst:freqhig;
s =freq;

freqpts =max(size(s));
for q = 1:freqpts,
    
        HL_ori = C_cyc_ori*((exp(sqrt(-1)*s(q))*E_cyc_ori-A_cyc_ori)\B_cyc_ori);
%         HL = C_cyc*((exp(sqrt(-1)*s(q))*E_cyc-A_cyc)\B_cyc);
        Hl = C_r*((exp(sqrt(-1)*s(q))*E_r-A_r)\B_r);
    
        nrmHL_ori(q)=norm(full(HL_ori));
%         nrmHL(q)=norm(full(HL));
        nrmHl(q)=norm(full(Hl));
    
        err(q)=norm(HL_ori-Hl);
           
end
% 
% figure(4); clf; % norms in log sclae
% loglog(s,nrmHL_ori,'r-o','linewidth',1);
% hold on
% loglog(s, nrmHl,'b:*','linewidth',1);
% hold off
% xlabel('Frequency (rad/sec)')
% ylabel('Magnitude')
% legend('Original System','Projected System','Reduced System')
% legend('Original System','Reduced System')

figure(5);  % norms in linear sclae
semilogx(s,nrmHL_ori,'r:o','linewidth',2);
hold on
semilogx(s, nrmHl,'b-*','linewidth',1);
hold off

xlabel('Frequency (rad/sec)')
ylabel('Magnitude')
legend('Original System','Reduced System')
%         
figure(6); 
semilogx(s,err,'b', s, ones(1,freqpts)*errb,'r--');
% semilogy(s,err,'b', s, errb,'r--');
xlabel('Frequency (rad/sec)')
ylabel('Absolute error')
legend('Error of System','Error Bound')

