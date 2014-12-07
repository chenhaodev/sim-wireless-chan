function CV(system_para);

n = 0;
system_para.t = log2(system_para.P0/system_para.P1);

while(n < system_para.num_of_test)
n = n + 1;

% step1 将alpha由高到低排序
system_para.alpha = system_para.alpha_d * system_para.Qd - system_para.alpha_f * system_para.Qf;
[system_para.alpha, IX] = sort( system_para.alpha, 'descend');

% step2 generate u(i), i = 1,..,M
PU_existance = randsrc(1,1,[-1,1;system_para.P0,system_para.P1]);

for i = 1:1:system_para.M
if(PU_existance == 1)
system_para.u(i) = randsrc(1,1,[1,-1;system_para.Pd(i),1 - system_para.Pd(i)]);
else
system_para.u(i) = randsrc(1,1,[1,-1;system_para.Pf(i),1 - system_para.Pf(i)]);
end
end
system_para.u = system_para.u(IX); % 将u按关键字alpha排序,靠前的u(i)对应的alpha较大

% step3
U = 0;

% sub_step1: get test statistics U
for i = 1:1:system_para.M
if(system_para.u(i) == 1)
U = U + log2( system_para.Qd(i) / system_para.Qf(i) );
else
U = U + log2( (1 - system_para.Qd(i)) / (1 - system_para.Qf(i)) );
end
end

% sub_step2: judgement, update PD and PF, N1 and N0
if( U >= system_para.t)
u_global = 1;
else
u_global = -1;
end

if(u_global == 1 && PU_existance == 1)
system_para.H1_num = system_para.H1_num + 1;
system_para.u1_num = system_para.u1_num + 1;
system_para.N1 = system_para.N1 + 1;

elseif(u_global == -1 && PU_existance == 1) % 漏警
system_para.H1_num = system_para.H1_num + 1;
system_para.N0 = system_para.N0 + 1;

elseif(u_global == 1 && PU_existance == -1) % 虚警
system_para.H0_num = system_para.H0_num + 1;
system_para.u0_num = system_para.u0_num + 1;
system_para.N1 = system_para.N1 + 1;

elseif(u_global == -1 && PU_existance == -1)
system_para.H0_num = system_para.H0_num + 1;
system_para.N0 = system_para.N0 + 1;
end

system_para.PD = system_para.u1_num / system_para.H1_num;
system_para.PF = system_para.u0_num / system_para.H0_num;


% sub_step3: update Qd and Qf
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

% sub_step4: update QD and QF
m0 = 0;
m1 = 0;
sigma0 = 0;
sigma1 = 0;

for i = 1:1:system_para.M
A = log2( system_para.Qd(i) / system_para.Qf(i) );
B = log2((1 - system_para.Qd(i)) / (1 - system_para.Qf(i)));
m1 = m1 + system_para.Qd(i) * A + (1 - system_para.Qd(i)) * B;
m0 = m0 + system_para.Qf(i) * A + (1 - system_para.Qf(i)) * B;
sigma1 = sigma1 + system_para.Qd(i) * (1 - system_para.Qd(i)) * ( A - B )^2;
sigma0 = sigma0 + system_para.Qf(i) * (1 - system_para.Qf(i)) * ( A - B )^2;
end

%    用Pd,Pf代替Qd,Qf，效果也不好
%
%      for i = 1:1:system_para.M
%          A = log2( system_para.Pd(i) / system_para.Pf(i) );
%          B = log2((1 - system_para.Pd(i)) / (1 - system_para.Pf(i)));
%          m1 = m1 + system_para.Pd(i) * A + (1 - system_para.Pd(i)) * B;
%          m0 = m0 + system_para.Pf(i) * A + (1 - system_para.Pf(i)) * B;
%          sigma1 = sigma1 + system_para.Pd(i) * (1 - system_para.Pd(i)) * ( A - B )^2;
%          sigma0 = sigma0 + system_para.Pf(i) * (1 - system_para.Pf(i)) * ( A - B )^2;
%      end

system_para.QD = 1 - normcdf( (system_para.t - m1) / sigma1 );
system_para.QF = 1 - normcdf( (system_para.t - m0) / sigma0 );

display( 'program is running....just finished a test' );
log.n=n;
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



end % end of while

display( 'running totally finished..............' );
plot(system_para.record);
title('M vs n');