clear all,close all,clc

addpath('./tools');
addpath('./data');

%%  generate Data

data=load('chua.mat');
dt = (data.s(2,1) - data.s(1,1));

data.CH1V = wiener2(data.CH1V,[30 1]);
data.CH2V = wiener2(data.CH2V,[30 1]);
N=50;

U = data.CH1V(mod(0:length(data.CH1V)-1,1+N)<1);
I = data.CH2V(mod(0:length(data.CH2V)-1,1+N)<1);

t = 1:size(U,1)*1.02e-7;
[U, dU, index] = estimatediff(U, t, 'solver', 1, []);
I = I(index);

index = 31:207;%because of wiener
U = U(index);
I = I(index);

dU = dU(index);

ans_sys_idx{1} = [];
ans_sys_idx{2} = [];
ans_sys_idx{3} = [];
for i = 1:size(U,1)
    if U(i,1)<-1.5&&U(i,1)>-6
        ans_sys_idx{1} = [ans_sys_idx{1} i ];
        state(i) = 1;
    end
    if U(i,1)<1.5&&U(i,1)>-1.5
        ans_sys_idx{2} = [ans_sys_idx{2} i ];
        state(i) = 2;
    end
    if U(i,1)<6&&U(i,1)>1.5
        ans_sys_idx{3} = [ans_sys_idx{3} i ];
        state(i) = 3;
    end
end


%%


A = [ones(size(U)) U I exp(U) U./I  cos(0.1*U).^2./(1+I.^2)  cos(U+I).^2];
% A = [ones(size(U)) U I ];

parameter.lambda = [0.05 0.01];    % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [0.012  0.044];
parameter.state = state;
parameter.Phi = A;
parameter.y = dU;
parameter.normalize_y = 1;
[result]=ihyde(parameter);



result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;

%%
% % Phi2 = [ ones(size(U)) U sin(I).*cos(U)  dU./(sin(U)+dU) dU./I  dU];
% 
% Phi2 = [ ones(size(U)) U ];
% 
% para_log.idx_sys = idx_sys;
% para_log.beta = 0.01;
% 
% para_log.y = dU;
% para_log.Phi2 = Phi2;
% 
% para_log.normalize = 1;
% 
% [syslogic,labelMat,data] = ihydelogic(para_log);
% 
% for i =1:size(para_log.idx_sys,2)
%     for j=1:size(para_log.idx_sys,2)
%         if length(syslogic{i,j})~=0
%             logicsys(i,j) = -syslogic{i,j}(1)/syslogic{i,j}(2);
%         end
%     end
%     
% end

%%
judge = zeros(size(state));
judge(idx_sys{1}) =3;
judge(idx_sys{2}) =1;
judge(idx_sys{3}) =2;
wrong_numbers = find((judge - state)~=0);
%%
% figure(1)
% axes1 = axes('Parent',figure(1));
% hold on
% color = {'r' ,[0.85 0.33 0.1],'b'};
% for i =[1,3,2]
%     input1 = zeros(size(U));
%     input1(ans_sys_idx{i},1) = U(ans_sys_idx{i},1);
%     
%     input1(input1==0)=nan;
%     
%     plot(input1(:,1),'Color',color{i},'LineWidth',2,'Marker','o','LineStyle','none');
%     legend('Subsystem_1','Subsystem_2','Subsystem_3')
%     
% end
% 
% 
% xlabel('Time(10^{-5}s)','FontWeight','bold');
% ylabel('y_1(V)','FontWeight','bold');
% box(axes1,'on');
% set(axes1,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
% legend(axes1,'show');
% %%
% figure(2)
% axes1 = axes('Parent',figure(2));
% hold on
% color = {'r' ,[0.85 0.33 0.1],'b'};
% for i =[1,3,2]
%     input1 = zeros(size(U));
%     input1(ans_sys_idx{i},1) = U(ans_sys_idx{i},1);
%     
%     input1(input1==0)=nan;
%     
%     plot(input1(:,1),'Color',color{i},'LineWidth',2,'Marker','o','LineStyle','none');
%     
%     
% end
% color = {'b' ,'r',[0.85 0.33 0.1]};
% for i =[2 1 3];
%     input1 = zeros(size(U));
%     input1(idx_sys{i},1) = U(idx_sys{i},1);
%     
%     input1(input1==0)=nan;
%     
%     plot(input1(:,1),'Color',color{i},'LineWidth',2,'Marker','+','LineStyle','none');
%     
%     
% end
% legend('Subsystem_1','Subsystem_2','Subsystem_3','Identified subsystem_1','Identified subsystem_2','Identified subsystem_3')
% xlabel('Time(10^{-5}s)','FontWeight','bold');
% ylabel('y_1(V)','FontWeight','bold');
% box(axes1,'on');
% set(axes1,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
% legend(axes1,'show');

