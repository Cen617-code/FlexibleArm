function kx = jiekx0new11(h1, h2, h3, ~, ~, ~, dh1, dh2, dh3, dt, q)
%JIEKX0NEW11 约束雅可比-速度项的数值求导
%   kx = jiekx0new11(...) 通过对 pqqd11 的有限差分，构造约束
%   与速度向量 q 的耦合项。

%% 提取速度分量
velocity = q(1:3);

%% 有限差分求导
baseValue = pqqd11(h1, h2, h3, velocity(1), velocity(2), velocity(3), dh1, dh2, dh3, dt, q);

perturbations = [dh1, 0, 0;
				 0, dh2, 0;
				 0, 0, dh3];

derivatives = zeros(6, 3);

for idx = 1:3
	perturbedH = [h1; h2; h3] + perturbations(:, idx);
	perturbedValue = pqqd11(perturbedH(1), perturbedH(2), perturbedH(3), ...
							velocity(1), velocity(2), velocity(3), ...
							dh1, dh2, dh3, dt, q);
	derivatives(:, idx) = (perturbedValue - baseValue) / perturbations(idx, idx);
end

%% 组装输出
kx = [derivatives, zeros(6, 6)] * q;
end