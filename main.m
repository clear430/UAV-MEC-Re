%Zhou F, Wu Y, Hu R Q, et al. 
%Computation Rate Maximization in UAV-Enabled Wireless-Powered Mobile-Edge Computing Systems[J]. 
%IEEE Journal on Selected Areas in Communications, 2018, 36(9): 1927-1941.

%Author: ZhuoYi Bai
%E-mail: oliverpai@163.com
%--------------------------�̶���������-----------------------------
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
M = 4;                                  %�����û���
Qm = zeros(M,2);                        %�����û�����
Wm = ones(M);                           %�����û�Ȩ�� 
N = 50;                                 %ʱ϶��
T = 2;                                  %��ʱ��2s
Height = 10;                            %���˻��߶�10m
C = 1000;                               %ÿbit��Ҫ1000����
Eta_0 = 0.8;                            %����ת��Ч��
B = 4*10^1;                             %����40MHz
Sigma = 10^-9;                          %��������10^-9W
Gamma_c = 10^-28;                       %������оƬ��Ч����ϵ��
Beta = 10^(-50/10);                     %��λ1m���ŵ���������-50dB
Epsilon = 10^-4;                        %�ݴ�
V_max = 20;                             %���˻��������ٶ�20m/s
P_0 = 0.1;                              %���˻����书��0.1W
vm = 1.1;                               %������������Ч��������ֵ����ӳ���⿪֧��

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
% ��ȡ·��
% [f_partial,P_partial,t_partial,Qu_partial,R_partial] = TwoStageAlgorithm(M,Wm,Qm);
% scatter(Qu_partial(:,1),Qu_partial(:,2),'g')
% Qu_semi = zeros(N,2);
% diamet = norm(Qm(M,:)-Qm(1,:))/2;
% theta = pi:-pi/(N-1):0;
% Qu_semi(:,1) = diamet*cos(theta)'+(Qm(M,1)-Qm(1,1))/2;
% Qu_semi(:,2) = diamet*sin(theta)'+(Qm(M,2)-Qm(1,2))/2;
% hold on;
% scatter(Qu_semi(:,1),Qu_semi(:,2),'r')

% ��������켣�������˻����书�ʶ������������
for i = 1:10
    P_0 = 0.1*i;
    [f_partial,P_partial,t_partial,Qu_partial,R_partial] = TwoStageAlgorithm(M,Wm,Qm);
    [f_semi,P_semi,t_semi,R_semi] = SemiCircle_R(M,Wm,Qm);
    Algorithm_R_p(i) = 10^6*sum(R_partial);
    Semi_R_p(i) = 10^6*sum(R_semi);
    display(['���书�ʣ�',num2str(P_0)]);
    display(['��������',num2str(Algorithm_R_p(i))]);
    display(['��������',num2str(Semi_R_p(i))]);
end
plot(0.1:0.1:1,Algorithm_R_p);
hold on;
plot(0.1:0.1:1,Semi_R_p);

%����ж��ģʽ




