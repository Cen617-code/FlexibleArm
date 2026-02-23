function beta = cosbeta1hanshu(h1, h2, h3, L0)
%COSBETA1HANSHU 计算柔性臂三角机构中的摆角 beta
%   beta = cosbeta1hanshu(h1, h2, h3, L0) 基于三个广义位移和底边长度
%   L0 计算关节平面内的摆角。返回值范围在 [-pi, pi]。

%% 预检查
if L0 <= 0
    error('cosbeta1hanshu:InvalidBaseLength', 'Base length L0 must be positive.');
end

sqrt3 = sqrt(3);

delta12 = h2 - h1;
delta13 = h3 - h1;

%% 定义几何向量
p1p2 = [sqrt3 * L0 / 2;
         3 * L0 / 2;
         delta12];

xita1 = [0;
          sqrt3 * delta12 * L0;
         -3 * sqrt3 * L0^2 / 2];

xita2 = [3 * delta13 * L0 / 2;
          sqrt3 * (delta12 * L0 - delta13 * L0 / 2);
         -3 * sqrt3 * L0^2 / 2];

%% 角度判定
crossVec = cross(xita1, xita2);
dotValue = dot(xita1, xita2);

normProduct = norm(xita1) * norm(xita2);
if normProduct < eps
    error('cosbeta1hanshu:DegenerateConfiguration', ...
          'The vectors defining beta are nearly linearly dependent.');
end

cosBeta = max(min(dotValue / normProduct, 1), -1);

if dot(p1p2, crossVec) >= 0
    beta = acos(cosBeta);
else
    beta = -acos(cosBeta);
end

beta = real(beta);
end
