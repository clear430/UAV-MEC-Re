%Zhou F, Wu Y, Hu R Q, et al. 
%Computation Rate Maximization in UAV-Enabled Wireless-Powered Mobile-Edge Computing Systems[J]. 
%IEEE Journal on Selected Areas in Communications, 2018, 36(9): 1927-1941.

%Author: ZhuoYi Bai
%E-mail: oliverpai@163.com
%--------------------------固定参数设置-----------------------------
global M;
global N;
global T;
global Height;
global C;
global Eta_0;
global B;
global Sigma_0;
global Gamma_c;
global Beta_0;
global Epsilon;
global Epsilon_1;
global Epsilon_2;
global V_max;
global P_0;
global vm;
global step_lambda;
global step_alpha;
M = 4;                                  %地面用户数
Qm = zeros(M,2);                        %地面用户坐标
Wm = ones(M);                           %地面用户权重 
N = 50;                                 %时隙数
T = 2;                                  %总时长2s
Height = 10;                            %无人机高度10m
C = 1000;                               %每bit需要1000周期
Eta_0 = 0.8;                            %能量转换效率
B = 40;                                 %带宽40MHz
Sigma_0 = 10^-9;                        %噪声功率10^-9W
Gamma_c = 10^-28;                       %处理器芯片有效电容系数
Beta_0 = 10^(-50/10);                   %单位1m处信道功率增益-50dB
Epsilon = 10^-4;                        %容错
Epsilon_1 = 10^-4;
Epsilon_2 = 10^-4;
V_max = 20;                             %无人机最大飞行速度20m/s
P_0 = 0.1;                              %无人机传输功率0.1W
vm = 1;                                 %总数据量与有效数据量比值（反映额外开支）
step_lambda = 0.1;
step_alpha = 0.1;

%-------------------------------初始化-----------------------------
Wm(1) = 0.1;
Wm(2) = 0.4;
Wm(3) = 0.3;
Wm(4) = 0.2;
Qm(1,:) = [0 0];                             
Qm(2,:) = [0 10];
Qm(3,:) = [10 10];
Qm(4,:) = [10 0];

%------------------------执行部分卸载模式算法-----------------------
[f_partial,P_partial,t_partial,Qu_partial] = TwoStageAlgorithm(Wm,Qm);

%二分卸载模式




