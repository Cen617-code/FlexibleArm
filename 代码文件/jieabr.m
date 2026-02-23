function [alpha, b, r] = jieabr(h1, h2, h3, L0)
%JIEABR 计算柔性臂的三角广义角度 (alpha, b, r)
%   [alpha, b, r] = jieabr(h1, h2, h3, L0) 返回根据广义位移 h1~h3
%   以及基座边长 L0 得到的角度参数。

%% 预计算量
if L0 <= 0
	error('jieabr:InvalidBaseLength', 'Base length L0 must be positive.');
end

delta12 = h2 - h1;
delta13 = h3 - h1;
delta23 = h3 - h2;

halfBase = 3 * L0 / 2;

%% 轴线旋转角
alpha = atan2(delta12, halfBase);

%% 节间摆角
b = cosbeta1hanshu(h1, h2, h3, L0);

%% 边长推导
L1 = sqrt(halfBase^2 + delta12^2) - L0 / 2;

%% 末端旋转角
L13Squared = 3 * L0^2 + delta13^2;
L23Squared = 3 * L0^2 + delta23^2;

[r, ~, ~] = gammafast(h1, h2, h3, L13Squared, L23Squared, L0, L1);
end