function constraintProduct = pqqd11(h1, h2, h3, ~, ~, ~, dh1, dh2, dh3, ~, qd)
%PQQD11 约束雅可比与状态向量的乘积 (数值求导实现)
%   constraintProduct = pqqd11(...) 通过有限差分计算位姿与角度对
%   广义坐标的偏导，进而与状态 qd 相乘得到约束项贡献。

%% 几何参数
geometry.h0 = 3;
geometry.L0 = 2;

%% 基础状态
h = [h1; h2; h3];
deltaH = [dh1; dh2; dh3];

poseCurrent = evaluatePose(h, geometry);
anglesCurrent = evaluateAngles(h, geometry);

%% 数值雅可比
positionJacobian = zeros(3);
angleJacobian = zeros(3);

for idx = 1:3
	step = deltaH(idx);
	if abs(step) < eps
		error('pqqd11:ZeroFiniteDifferenceStep', ...
			  'Finite difference step for coordinate %d must be non-zero.', idx);
	end

	perturbedH = h;
	perturbedH(idx) = perturbedH(idx) + step;

	posePerturbed = evaluatePose(perturbedH, geometry);
	anglesPerturbed = evaluateAngles(perturbedH, geometry);

	positionJacobian(:, idx) = (posePerturbed.position - poseCurrent.position) / step;
	angleJacobian(:, idx) = (anglesPerturbed - anglesCurrent) / step;
end

%% 组装雅可比矩阵并与状态向量相乘
constraintJacobian = [positionJacobian, -eye(3), zeros(3);
					  angleJacobian,    zeros(3), -eye(3)];

constraintProduct = constraintJacobian * qd;
end

%% ===================== 辅助函数 =====================

function pose = evaluatePose(h, geometry)
	transformation = jieA(h(1), h(2), h(3), geometry.h0, geometry.L0);
	originHomogeneous = transformation * [0; 0; 0; 1];

	pose.position = originHomogeneous(1:3);
	pose.rotation = transformation(1:3, 1:3);
end

function angles = evaluateAngles(h, geometry)
	[alpha, beta, gamma] = jieabr(h(1), h(2), h(3), geometry.L0);
	angles = [alpha; beta; gamma];
end