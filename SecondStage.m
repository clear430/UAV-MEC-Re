
%--------------------------����SecondStage-------------------------
%���ܣ�����P4Ŀ�꺯��
%���룺
%       P:                          �Ż����ʽ��
%       t:                          �Ż�ʱ����
%       q_last:                     ��һ�ε��������˻�����
%       q_var:                      �Ż����������ε������˻�����
%       W_m:                        �û�Ȩ������
%       Q_m:                        �û��������
%�����
%       p4_output:                  Ŀ�꺯�����
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
