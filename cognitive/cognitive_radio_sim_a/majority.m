function majority(system_para);


% 超过半数认为PU存在，fusion结果就是存在
%
% -------------------------------------------------------------------------
%  值得注意：k和threshold的取值对结果影响很大。
%  当 k = floor(system_para.M / 2) - 1，或ceil(system_para.M / 2) - 1，同时
%  threshold固定为0时，效果比较好
%
%  当 k = ceil(system_para.M / 2) - 1, threshold = 2 * k - system_para.M时，
%  运行结果不稳定 每次收敛的结果都不同
% -------------------------------------------------------------------------


n = 0;
system_para.t = 0;

system_para.Pd = ones(1, system_para.M) * 0.6; % 单个sensor的检测概率
system_para.Pf = ones(1, system_para.M) * 0.4; % 单个sensor的虚警概率

while(n < system_para.num_of_test)
n = n + 1;

k = ceil(system_para.M / 2); % k的取法是有讲究的，参见"distributed bayesian detection"

% step1 将alpha由高到低排序

% step2 generate u(i), i = 1,..,M
PU_existance = randsrc(1,1,[-1,1;system_para.P0,system_para.P1]);

for i = 1:1:system_para.M
if(PU_existance == 1)
system_para.u(i) = randsrc(1,1,[1,-1;system_para.Pd(i),1 - system_para.Pd(i)]);
else
system_para.u(i) = randsrc(1,1,[1,-1;system_para.Pf(i),1 - system_para.Pf(i)]);
end
end

% step3
% sub_step1 judgement and update PD, PF
U = 0;

for i = 1:1:system_para.M
U = U + system_para.u(i);
end

threshold = 2 * k - system_para.M;
if(U >= 0) % 0比threshold效果好
u_global = 1;
else
u_global = -1;
end

if(u_global == 1 && PU_existance == 1)
system_para.H1_num = system_para.H1_num + 1;
system_para.u1_num = system_para.u1_num + 1;
system_para.N1 = system_para.N1 + 1;

elseif(u_global == -1 && PU_existance == 1) %漏警
system_para.H1_num = system_para.H1_num + 1;
system_para.N0 = system_para.N0 + 1;

elseif(u_global == 1 && PU_existance == -1) %虚警
system_para.H0_num = system_para.H0_num + 1;
system_para.u0_num = system_para.u0_num + 1;
system_para.N1 = system_para.N1 + 1;

elseif(u_global == -1 && PU_existance == -1)
system_para.H0_num = system_para.H0_num + 1;
system_para.N0 = system_para.N0 + 1;
end

system_para.PD = system_para.u1_num / system_para.H1_num;
system_para.PF = system_para.u0_num / system_para.H0_num;

% sub_step2: update Qd and Qf

for i = 1:1:system_para.M
if( system_para.u(i) == 1 && u_global == 1 )
system_para.Qd(i) = ( system_para.Qd(i) * (system_para.N1 - 1) + 1 ) / system_para.N1;

elseif( system_para.u(i) == -1 && u_global == 1 )
system_para.Qd(i) = ( system_para.Qd(i) * (system_para.N1 - 1) )/ system_para.N1;

elseif( system_para.u(i) == 1 && u_global == -1 )
system_para.Qf(i) = ( system_para.Qf(i) * (system_para.N0 - 1) + 1 )/ system_para.N0;

elseif( system_para.u(i) == -1 && u_global == -1 )
system_para.Qf(i) = ( system_para.Qf(i) * (system_para.N0 - 1) )/ system_para.N0;
end
end

% sub_step4: update QD, QF
system_para.QD = 0;
system_para.QF = 0;
Qd_average = sum( system_para.Qd(1:system_para.M) ) / system_para.M;
Qf_average = sum( system_para.Qf(1:system_para.M) ) / system_para.M;

for i = k:1:system_para.M
A = factorial(system_para.M) / factorial(i) / factorial(system_para.M - i);
B = Qd_average^i * ( 1 - Qd_average )^(system_para.M - i);
C = Qf_average^i * ( 1 - Qf_average )^(system_para.M - i);
system_para.QD = system_para.QD + A * B;
system_para.QF = system_para.QF + A * C;
end


display( 'program is running....just finished a test' );
log.n=n;
log.PD_t=system_para.PD_t;
log.PF_t=system_para.PF_t;
log.Qd_average=Qd_average;
log.Qf_average=Qf_average;
log.M=system_para.M;
log.QD=system_para.QD;
log.QF=system_para.QF;
log.PD=system_para.PD;
log.PF=system_para.PF;
log
system_para.record(n) = system_para.M;


% sub_step5: update M
if(system_para.QD > system_para.PD_t + system_para.e1 && system_para.QF < system_para.PF_t - system_para.e1)
system_para.M = system_para.M - system_para.step;
elseif(system_para.QD < system_para.PD_t - system_para.e2 || system_para.QF > system_para.PF_t + system_para.e2)
system_para.M = system_para.M + system_para.step;
end

if(system_para.M > system_para.N)
system_para.M = system_para.N;
elseif(system_para.M <= 1)
system_para.M = 1;
end

end

display( 'running totally finished..............' );
plot(system_para.record);
title('M vs n');