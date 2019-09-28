function rho = modified_rho_fun(T,ef,Ib2)
% 单位均为国际单位
me = 9.10956e-31; % 电子质量kg
qe = 1.602e-19; % 元电荷C
mu0 = 4*pi*1e-7; % 真空磁导率
N = 525; % 螺线管匝数
D = 0.065; % 螺线管平均直径
L = 0.145; % 螺线管平均长度
d = 8.4e03; % 圆柱面极的直径
r = d/2; % 半径
K = 1.0709974e-18; % 能到达的最小动能ek=K*Ib2
ek = Ib2; %
B = mu0*N*Ib2*K/sqrt(L^2+D^2);
R = sqrt(2*ek*K*me)/qe/B;
if r/2 < R && R < r
    alpha = pi - asin(d/2/R);
elseif R > r
    alpha = asin(d/2/R);
else
    alpha = inf; % 无法到达
end
g = @(T,ef,x) exp((x-ef)*K/k/T)*K./( k*T* (exp((x-ef)*K/k/T)+1).^2 ); % 费米分布的密度函数
rho = 1./alpha * g(T,ef,ek);