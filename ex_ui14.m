clear
clc
close all

addpath('./tools');
addpath('./data');
%%  load Data
data = load('TwoMissLine.mat');
data1 = data.allData{1};
data2 = data.allData{2};
u1 = data1{1}';
i1 = data1{2}';
y1 = data1{3}';

u2 = data2{1}';
i2 = data2{2}';
y2 = data2{3}';

index = 1:10;
u2 = u2(index,:);
i2 = i2(index,:);


%% library
A1 = [real(u1) -imag(u1);imag(u1) real(u1) ];
A2 = [real(u2) -imag(u2);imag(u2) real(u2) ];



A = [A1;A2];
ii1 = [real(i1);imag(i1)];
ii2 = [real(i2);imag(i2)];

i = [ii1;ii2];
%% identify subsystem
choose = 6;
i = i(:,choose);
parameter.lambda = [0.001 0.5e-5];  % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [0.007  0.05];



parameter.Phi = A;
parameter.y = i;
parameter.normalize_y = 1;
[result]=ihyde(parameter);


result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;

%% identify logic
index = 1:size(A,1);
index = index';
Phi2 = [ones(length(index),1) index i i.^2 index./(index+10*i)];


para_log.idx_sys = idx_sys;
para_log.beta = 0.05;

para_log.y = i;
para_log.Phi2 = Phi2;

para_log.normalize = 1;
[syslogic,labelMat,data] = ihydelogic(para_log);