
%-----------------------SemiCircle_R----------------------
%功能：半圆轨迹部分卸载模式下的优化算法
%输入：
%       M:                               用户数量
%       W_m:                             用户权重向量
%       Q_m:                             用户位置向量
%输出：    
%       f_opt：                          各时隙最优CPU周期数
%       P_opt：                          各时隙最优传输功率
%       t_opt：                          各时隙最优卸载时间
%       qu_opt：                         各时隙无人机最优位置
%       R_opt：                          最优总计算量
function[f_opt,P_opt,t_opt,R_opt] = SemiCircle_R(M,W_m,Q_m)
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
global P_0;
global vm;
    f_opt = zeros(M,N);
    P_opt = zeros(M,N);
    t_opt = zeros(M,N);
    Qu_opt = zeros(N,2);
    diamet = norm(Q_m(M,:)-Q_m(1,:))/2;
    theta = pi:-pi/(N-1):0;
    Qu_opt(:,1) = diamet*cos(theta)'+(Q_m(M,1)-Q_m(1,1))/2;
    Qu_opt(:,2) = diamet*sin(theta)'+(Q_m(M,2)-Q_m(1,2))/2;
%    Qu_opt(:,1) = Q_m(1,1):(Q_m(M,1)-Q_m(1,1))/(N-1):Q_m(M,1);
    H_mn = zeros(M,N);
    R_sum_last = 0;
    while(1)
        cvx_begin
            cvx_expert true
            variable f_temp(M,N)
            variable Z_temp(M,N)
            variable t_temp(M,N)
            expression p2_obj(M,N)
            expression c5_left(M,N)
            c5_right = zeros(M,N);
            for user = 1:M
                for slot = 1:N
                    H_mn(user,slot) = Beta/(Height^2+norm(Qu_opt(slot,:)-Q_m(user,:))^2);
                    p2_obj(user,slot) = W_m(user)*...
                        (T*f_temp(user,slot)/N/C+...
                        B*T*(-rel_entr(t_temp(user,slot),t_temp(user,slot)+H_mn(user,slot)/Sigma*Z_temp(user,slot)))/(log(2)*vm*N));
                end
            end
            c5_right(:,1) = Eta_0*T/N*P_0*H_mn(:,1);
            for slot = 2:N
                c5_right(:,slot) = c5_right(:,slot-1)+Eta_0*T/N*P_0*H_mn(:,slot);
            end
            maximize sum(sum(p2_obj))
            subject to
                c5_left(:,1) = Gamma_c*pow_p(f_temp(:,1),3)+Z_temp(:,1);
                for slot = 2:N
                    c5_left(:,slot) = c5_left(:,slot-1)+Gamma_c*pow_p(f_temp(:,slot),3)+Z_temp(:,slot);
                end
                %c5_left = Gamma_c*pow_p(f_temp,3)+Z_temp;
                for slot = 1:N
                    for user = 1:M
                        f_temp(user,slot) >= 0;
                        Z_temp(user,slot) >= 0;
                        t_temp(user,slot) >= 0;
                          %每个时隙累积能量限制
                        %c5_left(user,slot)<= c5_right(user,slot);%Eta_0*T/N*P_0*H_mn(user,slot);
                        Gamma_c*pow_p(f_temp(user,slot),3)+Z_temp(user,slot) <= Eta_0*T/N*P_0*H_mn(user,slot);
                    end
                    sum(t_temp(:,slot)) <= 1;
                end
                %总能量限制
%                 for user = 1:M
%                     sum(c5_left(user,:)) <= Eta_0*T/N*P_0*sum(H_mn(user,:));
%                 end
        cvx_end
        f_opt = f_temp;
        P_opt = Z_temp./t_temp;
        t_opt = t_temp;
%计算总数据量
        R_sum = cvx_optval;
        display(['总数据量',num2str(R_sum)]);
        if (R_sum-R_sum_last) <= Epsilon
            R_opt = R_sum;
            break;
        end
        R_sum_last = R_sum;
    end
end