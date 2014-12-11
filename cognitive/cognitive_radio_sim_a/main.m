clear all;
close all;
randn( 'state', 0 );

option = 2;
% 0: chair-varshney rule;
% 1: chair-varshney rule with PU alerting;
% 2: majority rule, sensors have same Pd(i)
% 3: majority rule with PU alerting
%
% PU alerting是指，当PU存在，但是fusion center没有检测到，PU会发一个警报
% 给fusion center，这样fusion center 就确知有一次漏警


system_para.num_of_test = 1000;
system_para.record = zeros(1, system_para.num_of_test);

system_para.N = 100;
system_para.M = system_para.N;

system_para.step = 1;
system_para.H1_num = 1; %PU出现的次数
system_para.H0_num = 1; %PU没出现的次数
system_para.u1_num = 1; %PU出现时，fusion center成功检测到的次数
system_para.u0_num = 0; %PU没出现时，fusion center误检测到的次数

system_para.ui1_num = zeros(1,system_para.M); % fusion center认为PU出现时， sensor i也认为出现的次数
system_para.ui0_num = zeros(1,system_para.M); % fusion center认为PU没出现时， sensor i却认为出现的次数

system_para.N1 = 1; %fusion center认为PU出现的次数
system_para.N0 = 1; %fusion center认为PU没出现的次数

system_para.P0 = 0.8;
system_para.P1 = 0.2;

system_para.PD_t = 0.90;
system_para.PF_t = 0.15;

system_para.Pd = randint(1,system_para.M,[500,800]);
system_para.Pd = system_para.Pd/1000;
system_para.Pf = randint(1,system_para.M,[200,500]);
system_para.Pf = system_para.Pf/1000;

system_para.Qd_init = 0.7; % 对于majority rule，初始值的选取会影响收敛速度，但不影响收敛值
system_para.Qf_init = 0.3;

system_para.Qd = ones(1,system_para.M);
system_para.Qd = system_para.Qd * system_para.Qd_init;
system_para.Qf = ones(1,system_para.M);
system_para.Qf = system_para.Qf * system_para.Qf_init;

system_para.QD = 0;
system_para.QF = 0;
system_para.PD = 0;
system_para.PF = 0;

system_para.alpha_d = 1;
system_para.alpha_f = 0;

system_para.u = zeros(1,system_para.M);

system_para.e1 = 0.001;
system_para.e2 = 0.03;


if(option == 0)
CV(system_para);
elseif(option == 1)
CV_with_alerting(system_para);
elseif(option == 2)
majority(system_para);
elseif(option == 3)
majority_with_alerting(system_para);
end