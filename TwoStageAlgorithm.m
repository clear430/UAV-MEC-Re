
%-----------------------函数TwoStageAlgorithm----------------------
%功能：部分卸载模式下的优化算法
%输入：
%       W_m:                             用户权重向量
%       Q_m:                             用户位置向量
%输出：    
%       f_opt：                          各时隙最优CPU周期数
%       P_opt：                          各时隙最优传输功率
%       t_opt：                          各时隙最优卸载时间
%       qu_opt：                         各时隙无人机最优位置
function[f_opt,P_opt,t_opt,Qu_opt] = TwoStageAlgorithm(W_m,Q_m)
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
global V_max;
global P_0;
global vm;
global step_lambda;
global step_alpha;
    f_opt = zeros(M,N);
    P_opt = zeros(M,N);
    t_opt = 0.25*ones(M,N);
    Qu_opt = zeros(N,2);
    lambda = abs(randn(M,N));
    alpha = abs(randn(N));
    H_mn = zeros(M,N);
    for slot=1:N
        R_sum_last = 0;
        while(1)
            for user = 1:M
                H_mn(user,slot) = Beta_0/(Height^2+norm(Qu_opt(slot,:)-Q_m(user,:))^2);
                f_opt(user,slot)= sqrt(W_m(user)/(3*C*Gamma_c*sum(lambda(user,slot:N))));
                P_opt(user,slot)= max([0,W_m(user)*B/vm/log(2*sum(lambda(user,slot:N)))-Sigma_0^2/H_mn(user,slot)]);
                syms t_var;
                t_opt(user,slot)= solve(log2(1+H_mn(user,slot)*P_opt(user,slot)/Sigma_0^2)-...
                    H_mn(user,slot)*t_var*P_opt(user,slot)/log(2*(Sigma_0^2*t_var+H_mn(user,slot)*t_var*P_opt(user,slot)))...
                    -vm*N*alpha(slot)/B/T==0,t_var)
                delta_lambda = Eta_0*T/N*sum(P_0*H_mn(user,1:slot))-T/N*sum(Gamma_c*f_opt(user,1:slot)^3+...
                    t_opt(user,1:slot).*P_opt(user,1:slot));
                lambda(user,slot)=lambda(user,slot)-step_lambda*delta_lambda;
            end
            delta_alpha = 1-sum(t_opt(:,slot));
            alpha(slot)=alpha(slot)-step_alpha*delta_alpha;
            while(1)
                cvx_begin quiet
                    variable Qu_temp(2)
                    maximize(SecondStage(P_opt,t_opt,Qu_opt(slot,:),Qu_temp,W_m,Q_m))
                    subject to
                        norm(Qu_temp-Qu_opt(slot,:))^2 <= V_max*T/N;
                        Qu_opt(1,:) == Q_m(0,:);
                        Qu_opt(N,:) == Q_m(M,:);
                        for user = 1:M
                            sum(Gamma_c*f_opt(user,:).^3+t_opt(user,:).*P_opt(user,:))...
                                <= Eta_0*P_0*Beta_0*sum((Height^2+2*norm(Qu_opt(slot,:)-Q_m(user,:))^2-norm(Qu_temp-Q_m(user,:))^2)...
                                ./(Height^2+norm(Qu_opt(slot,:)-Q_m(user,:))^2));
                        end
                cvx_end
                if norm(Qu_opt(slot,:)-Qu_temp)^2 <= Epsilon
                    break
                end
                Qu_opt(slot,:) = Qu_temp;
            end
            R_sum = 0;
            for user = 1:M
                R_sum = R_sum+W_m(user)*...
                    (T*f_opt(user,slot)/N/C+...
                    B*T*t_opt(user,slot)/vm/N*log2(1+H_mn(user,slot)*P_opt(user,slot)/Sigma_0^2));
            end
            if (R_sum-R_sum_last) <= Epsilon_1
                break
            end
            R_sum_last = R_sum;
        end
    end
end
