
clear
clc
close all

addpath('./tools');
addpath('./data');
addpath('./savemat');


%% design


load('current.mat')

%% make dic

index = 1:300:length(current);
current = current(index);
state = state(index);

y = current(2:end);
state = state(2:end);
x = current(1:end-1);
polyorder = 2;
usesine = 0;
memory =5;
basis_function.work = 'off';
dic= library(x,polyorder,usesine,memory,basis_function);

A = dic(memory+2:end,:);
y = y(memory+2:end,:);

state = state(memory+2:end,:);

%% identify subsystem


parameter.lambda = [0.3 0.00001];   % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [1e-4  0.026];

parameter.Phi = A;
parameter.y = y;
parameter.normalize_y = 1;
[result]=ihyde(parameter);

% [check]=ihyde(dy(index),A(index,:),maxs, epsilon,z_lambda,MAXITER,normalize_dy);



result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];

result.mode = 'knn';
result.half_win_size = 5;
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;



%%
% Phi2 = [ones(size(x)) x ];
t = 1:length(y);

Phi2 = [ones(size(t')) t' ];

para_log.idx_sys = idx_sys;
para_log.beta = 0.1;

para_log.y = y;
para_log.Phi2 = Phi2;

para_log.normalize = 1;

[syslogic,labelMat,data] = ihydelogic(para_log);


%%
close all
figure
hold on
index1 = 1:length(state);
y11 = y(index1);
judge11 = state(index1);
x11 = x(index1);
A1 = A(index1,:);

ansy1 = 1000000*ones(size(y11(:,1)));
ansy1(find(judge11==0)) = y11(find(judge11==0));
ansy1(ansy1==1000000)=nan;


ansy11 = 1000000*ones(size(y11(:,1)));
ansy11(idx_sys{1},:) = A1(idx_sys{1},:)*sys(:,1);
ansy11(ansy11==1000000)=nan;

err1 = ansy1 - ansy11;

ansy2 = 1000000*ones(size(y11(:,1)));
ansy2(find(judge11==1)) = y11(find(judge11==1));
ansy2(ansy2==1000000)=nan;


ansy22 = 1000000*ones(size(y11(:,1)));
ansy22(idx_sys{2},:) = A1(idx_sys{2},:)*sys(:,2);
ansy22(ansy22==1000000)=nan;

err2 = ansy2 - ansy22;

plot(ansy2,'r','Linewidth',5)
plot(ansy22,'ro','MarkerSize',10)
plot(ansy1,'b','Linewidth',5);
plot(ansy11,'bo','MarkerSize',10)
legend('normal','subsystem1','fault','subsystem2')


