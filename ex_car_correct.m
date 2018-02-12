clear all,close all,clc

addpath('./tools');
addpath('./data');

%%  load Data



basis_function.work='off';
data=load('normal_car.mat');  %%  load Data
index = 1000:1400;
% index = 1:1934;
flag = data.flag(index);   % 1:straghtway   0:curve
dy = data.dy(index);
v = data.v(index);


%% library
polyorder = 1;
memory = 1;
usesine  = 0;
A= library(v,polyorder,usesine,memory,basis_function);
% v_k1 = v(memory+1:end-1,:);
% v_k2 = v(memory:end-2,:);
% v_k3 = v(memory-1:end-3,:);
v = v(memory+2:end,:);

A = A(memory+2:end,:);
dy = dy(memory+2:end);
flag = flag(memory+2:end);

%% identify subsystem

parameter.lambda = [0.1 1e-5];   % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [100  8];
% parameter.state = flag;
parameter.Phi = A;
parameter.y = dy;
parameter.normalize_y = 1;
[result]=ihyde(parameter);



result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;


%% identify logic
% Phi2 = [ones(size(flag)) flag 1./v sin(v) cos(v) v.^2  v_k1./v_k2 v_k3.^2 ];
Phi2 = [ones(size(flag)) flag ];

para_log.idx_sys = idx_sys;
para_log.beta = .5;
para_log.y = dy;
para_log.Phi2 = Phi2;
para_log.normalize = 1;

[syslogic,labelMat,data] = ihydelogic(para_log);

%% compare with the right answer
judge = 7*ones(size(flag));
judge(idx_sys{1}) = 0;
judge(idx_sys{2}) = 1;
wrong_position = find((judge-flag)~=0);
ans_sys_idx{1} = find(flag==1);
ans_sys_idx{2} = find(flag==0);
%%
close all
figure(1)
axes1 = axes('Parent',figure(1));
hold on
color = {'r' ,'b'};
for i =1:2
    input1 = zeros(size(v));
    input1(ans_sys_idx{i},1) = v(ans_sys_idx{i},1);
    
    input1(input1==0)=nan;
    
    plot(input1(:,1),'Color',color{i},'LineWidth',3);
    
    
end
legend('Subsystem_1','Subsystem_2')

xlabel('Time(4ms)','FontWeight','bold');
ylabel('Velocity(cm/s)','FontWeight','bold');
box(axes1,'on');
set(axes1,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
legend(axes1,'show');
%%
figure(2)
axes1 = axes('Parent',figure(2));
hold on
color = {'r' ,'b'};
for i =1:2
    input1 = zeros(size(dy));
    input1(ans_sys_idx{i},1) = dy(ans_sys_idx{i},1);
    
    input1(input1==0)=nan;
    
    plot(input1(:,1),'Color',color{i},'LineWidth',3);
    
    
end
legend('Subsystem_1','Subsystem_2')

xlabel('Time(4ms)','FontWeight','bold');
ylabel('du','FontWeight','bold');
box(axes1,'on');
set(axes1,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
legend(axes1,'show');


