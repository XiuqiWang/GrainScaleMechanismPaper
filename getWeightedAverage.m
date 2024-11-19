function weighted_mean=getWeightedAverage(A)
[PDF, x_values] = ksdensity(A);
weighted_mean = sum(x_values .* PDF) / sum(PDF);
end

function [mu, sigma, mean_value, std_value]=getFitData(A)
% 使用 fitdist 函数拟合对数正态分布
pd = fitdist(A, 'Lognormal');

% 获取拟合分布的参数
mu = pd.mu; % 对数正态分布的均值（在对数空间）
sigma = pd.sigma; % 对数正态分布的标准差（在对数空间）

% 计算原空间中的均值和标准差
mean_value = exp(mu + (sigma^2) / 2);
std_value = sqrt((exp(sigma^2) - 1) * exp(2*mu + sigma^2));
end

% for i=0:4
% for j=1:3
% eval(['[mu(',num2str(j),',',num2str(i),'+1), sigma(',num2str(j),',',num2str(i),'+1), mean_value(',num2str(j),',',num2str(i),'+1), std_value(',num2str(j),',',num2str(i),'+1)]=getFitData(V_E',num2str(i),num2str(j),');']);
% end
% end