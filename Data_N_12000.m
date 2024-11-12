
 

%%% This is an index 2 problem%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------
%stable data:::::::::::::::5 diagonal
clear all; close all; clc; 
%-----------------------------------------------
n= 5000;   %order of the system(or 500,1000);
l= 2000;   %algebraic part
nin= 3;   %(number of input)
nout= 2;  % number of output 
i=3;      %period
n_p= 2*n - l; % size of projectd system

I=speye(n);
M= I+spdiags(i*ones(n,1),i,n,n)+spdiags(-i*i*ones(n,1),i,n,n)+spdiags(i*ones(n,1),-i*i,n,n)+spdiags(-i*ones(n,1),i*i,n,n);
K_uu= spdiags(i*ones(n,1),-i*i,n,n)+spdiags(-i*ones(n,1),-i*i,n,n)+spdiags(-i*ones(n,1),i*i,n,n)+spdiags(-i*i*ones(n,1),i,n,n)+spdiags(i*i*ones(n,1),-i,n,n);
% D=cell(1,K);
mu= 0.5+.1*i;
nu= 0.8+.1*i;%here
D= mu*M+nu*K_uu;
den= 0.02;
K_pp= i*spdiags(i*ones(l,1),-i*i,l,l)+spdiags(-i*ones(l,1),i*i,l,l)+spdiags(i*i*ones(l,1),-i*i,l,l)+spdiags(-i*i*ones(n,1),i,l,l)+spdiags(i*ones(n,1),-i*i,l,l);
%mu=0.005;nu=1;
%K_up=spdiags(ones(n,1),10,n,n)
s1=rng
K_up= sprand(n,l,0.1*den);
rng(s1)
%%%%%%%%%%%%%Transformation to first order system%%%%%%%%%%%%%%%%%%%%%%%%%
E1 = 2*eye(2*n);
% E1=[I spalloc(n,n,0);spalloc(n,n,0) I];
% E1=[M spalloc(n,n,0);spalloc(n,n,0) spalloc(n,n,0)];
% E1=[M M';M' M];
A1=0.05*[M K_uu';K_uu D];
%J1=J10;
s2=rng
A2=0.05*[sprand(n,l,0.1);K_up];
rng(s2)
% A3{i}=[-K_up' sprand(l,n,0.1) ];
A3=A2';
% J4=-K_pp; 
A4=[spalloc(l,l,0)];


% B1=[zeros(n,1);ones(n,1)];
% B2=ones(l,1);
% C1=[ones(1,n) zeros(1,n)];
% C2=B2';
B1=[spalloc(n,nin,0);spdiags(ones(n,1),0,n,nin)];
B2=spdiags(zeros(l,1),0,l,nin);
C1=[spdiags(ones(n,1),0,nout,n) spalloc(nout,n,0)];
C2=spdiags(zeros(l,1),0,nout,l);
E= [E1 spalloc(size(A2,1),size(A2,2),0);spalloc(size(A3,1),size(A3,2),0) spalloc(size(A4,1),size(A4,2),0)];

A= [A1 A2;A3 A4];
B= [B1; B2];
C= [C1 C2];
figure(1)
spy(A);
figure(2)
spy(E);
 

clear M K_uu D K_pp K_up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 7 diagonal

%% creating nonsingular system matrices

% Pl = eye(size(E1,1)) - (A2 * inv(A2' * inv(E1) *A2) * A2' * inv(E1));
Pl=eye(size(E1,1)) - A2*((A2'*(E1\A2))\(A2'/E1));
%   
 E_m = Pl*E1*Pl';
% 
 A_m = Pl*A1*Pl';
% 
 B_m = Pl*B1;
% 
 C_m = C1*Pl';
% 



