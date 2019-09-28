close all
K = 1.0709974e-18;
k = 1.38e-23;
e = 1.602e-19;
T = 1500;
data = importdata('fit.csv');
ib2 = data(:,1);
ip = data(:,2);
delta_ib2 = data(2:end,1);
delta_ip = data(2:end,3);
%% 用自定义函数拟合费米分布
g_e  = @(T,ef,x) 1./(exp((x-ef)*K/k/T)+1); % 拟合函数
dg_e = @(T,ef,x) exp((x-ef)*K/k/T)*K./( k*T* (exp((x-ef)*K/k/T)+1).^2 );
e_k = ib2; % 横坐标, 粒子能量
n_e = ip; % 纵坐标，粒子分布
results = fit(e_k,n_e,g_e,'StartPoint',[1500,0.2]);
delta_n_e = delta_ip/0.04; % 作差除以间距才是导数
results2 = fit(delta_ib2,delta_n_e,dg_e,'StartPoint',[3500,0.2]);
figure
hold on
plot(results,'r-');
plot(e_k,n_e,'bo--');
legend('拟合结果','实验结果');
title('分布函数拟合')
saveas(gca,'./dist_fun_fit.png');
hold off
figure
hold on
plot(results2,'r-');
plot(delta_ib2,delta_n_e,'bo--');
legend('拟合结果','实验结果');
title('密度函数拟合');
saveas(gca,'./rho_fun_fit.png');


