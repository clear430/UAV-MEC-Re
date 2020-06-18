%Zhou F, Wu Y, Hu R Q, et al. 
%Computation Rate Maximization in UAV-Enabled Wireless-Powered Mobile-Edge Computing Systems[J]. 
%IEEE Journal on Selected Areas in Communications, 2018, 36(9): 1927-1941.

%Author: ZhuoYi Bai
%E-mail: oliverpai@163.com
%--------------------------固定参数设置-----------------------------
global N;
global T;
global Height;
global C;
global Eta_0;
global B;
global Sigma;
global Gamma_c;
global Beta;
global Epsilon;
global V_max;
global P_0;
global vm;
M = 4;                                  %地面用户数
Qm = zeros(M,2);                        %地面用户坐标
Wm = ones(M);                           %地面用户权重 
N = 50;                                 %时隙数
T = 2;                                  %总时长2s
Height = 10;                            %无人机高度10m
C = 1000;                               %每bit需要1000周期
Eta_0 = 0.8;                            %能量转换效率
B = 4*10^1;                             %带宽40MHz
Sigma = 10^-9;                          %噪声功率10^-9W
Gamma_c = 10^-28;                       %处理器芯片有效电容系数
Beta = 10^(-50/10);                     %单位1m处信道功率增益-50dB
Epsilon = 10^-4;                        %容错
V_max = 20;                             %无人机最大飞行速度20m/s
P_0 = 0.1;                              %无人机传输功率0.1W
vm = 1.1;                               %总数据量与有效数据量比值（反映额外开支）

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
% 获取路径
% [f_partial,P_partial,t_partial,Qu_partial,R_partial] = TwoStageAlgorithm(M,Wm,Qm);
% scatter(Qu_partial(:,1),Qu_partial(:,2),'g')
% Qu_semi = zeros(N,2);
% diamet = norm(Qm(M,:)-Qm(1,:))/2;
% theta = pi:-pi/(N-1):0;
% Qu_semi(:,1) = diamet*cos(theta)'+(Qm(M,1)-Qm(1,1))/2;
% Qu_semi(:,2) = diamet*sin(theta)'+(Qm(M,2)-Qm(1,2))/2;
% hold on;
% scatter(Qu_semi(:,1),Qu_semi(:,2),'r')

% 用来计算轨迹下随无人机发射功率而变的总数据量
for i = 1:10
    P_0 = 0.1*i;
    [f_partial,P_partial,t_partial,Qu_partial,R_partial] = TwoStageAlgorithm(M,Wm,Qm);
    [f_semi,P_semi,t_semi,R_semi] = SemiCircle_R(M,Wm,Qm);
    Algorithm_R_p(i) = 10^6*sum(R_partial);
    Semi_R_p(i) = 10^6*sum(R_semi);
    display(['发射功率：',num2str(P_0)]);
    display(['数据量：',num2str(Algorithm_R_p(i))]);
    display(['数据量：',num2str(Semi_R_p(i))]);
end
plot(0.1:0.1:1,Algorithm_R_p);
hold on;
plot(0.1:0.1:1,Semi_R_p);

%二分卸载模式




