%Zhou F, Wu Y, Hu R Q, et al. 
%Computation Rate Maximization in UAV-Enabled Wireless-Powered Mobile-Edge Computing Systems[J]. 
%IEEE Journal on Selected Areas in Communications, 2018, 36(9): 1927-1941.

%Author: ZhuoYi Bai
%E-mail: oliverpai@163.com
%--------------------------�̶���������-----------------------------
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
M = 4;                                  %�����û���
Qm = zeros(M,2);                        %�����û�����
Wm = ones(M);                           %�����û�Ȩ�� 
N = 50;                                 %ʱ϶��
T = 2;                                  %��ʱ��2s
Height = 10;                            %���˻��߶�10m
C = 1000;                               %ÿbit��Ҫ1000����
Eta_0 = 0.8;                            %����ת��Ч��
B = 40;                                 %����40MHz
Sigma_0 = 10^-9;                        %��������10^-9W
Gamma_c = 10^-28;                       %������оƬ��Ч����ϵ��
Beta_0 = 10^(-50/10);                   %��λ1m���ŵ���������-50dB
Epsilon = 10^-4;                        %�ݴ�
Epsilon_1 = 10^-4;
Epsilon_2 = 10^-4;
V_max = 20;                             %���˻��������ٶ�20m/s
P_0 = 0.1;                              %���˻����书��0.1W
vm = 1;                                 %������������Ч��������ֵ����ӳ���⿪֧��
step_lambda = 0.1;
step_alpha = 0.1;

%-------------------------------��ʼ��-----------------------------
Wm(1) = 0.1;
Wm(2) = 0.4;
Wm(3) = 0.3;
Wm(4) = 0.2;
Qm(1,:) = [0 0];                             
Qm(2,:) = [0 10];
Qm(3,:) = [10 10];
Qm(4,:) = [10 0];

%------------------------ִ�в���ж��ģʽ�㷨-----------------------
[f_partial,P_partial,t_partial,Qu_partial] = TwoStageAlgorithm(Wm,Qm);

%����ж��ģʽ




