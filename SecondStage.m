
%--------------------------函数SecondStage-------------------------
%功能：论文P4目标函数
%输入：
%       P:                          优化功率结果
%       t:                          优化时间结果
%       q_last:                     上一次迭代的无人机坐标
%       q_var:                      优化变量，本次迭代无人机坐标
%       W_m:                        用户权重向量
%       Q_m:                        用户坐标矩阵
%输出：
%       p4_output:                  目标函数结果
function[p4_output] = SecondStage(P,t,q_last,q_var,W_m,Q_m)
global M;
global N;
global T;
global Height;
global B;
global Sigma_0;
global Beta_0;
global vm;
    p4_output = 0;
    for user = 1:M
        p4_output = p4_output+W_m(user)*B*T*sum(t(user,:).*...
            (log2(1+Beta_0*P(user,:)/Sigma_0^2./(Height^2+norm(q_last-Q_m(user,:))^2))...
            -Beta_0*P(user,:)*log2(exp(1))./(Sigma_0^2*Height^2+Beta_0*P(user,:)+Sigma_0^2*norm(q_last)^2)./(Height^2+norm(q_last)^2)...
            .*(norm(q_var)^2-norm(q_last)^2)))...
            /vm/N;
    end
end
