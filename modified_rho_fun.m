function rho = modified_rho_fun(T,ef,Ib2)
% ��λ��Ϊ���ʵ�λ
me = 9.10956e-31; % ��������kg
qe = 1.602e-19; % Ԫ���C
mu0 = 4*pi*1e-7; % ��մŵ���
N = 525; % ���߹�����
D = 0.065; % ���߹�ƽ��ֱ��
L = 0.145; % ���߹�ƽ������
d = 8.4e03; % Բ���漫��ֱ��
r = d/2; % �뾶
K = 1.0709974e-18; % �ܵ������С����ek=K*Ib2
ek = Ib2; %
B = mu0*N*Ib2*K/sqrt(L^2+D^2);
R = sqrt(2*ek*K*me)/qe/B;
if r/2 < R && R < r
    alpha = pi - asin(d/2/R);
elseif R > r
    alpha = asin(d/2/R);
else
    alpha = inf; % �޷�����
end
g = @(T,ef,x) exp((x-ef)*K/k/T)*K./( k*T* (exp((x-ef)*K/k/T)+1).^2 ); % ���׷ֲ����ܶȺ���
rho = 1./alpha * g(T,ef,ek);