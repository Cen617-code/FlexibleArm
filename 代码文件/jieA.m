function T = jieA(h1, h2, h3, h0, L0)
%JIEA 构造柔性臂的齐次变换矩阵
%   T = jieA(h1, h2, h3, h0, L0) 返回基于广义坐标 h1~h3 的平台位姿矩阵。

%% 基础量计算
delta12 = h2 - h1;
delta13 = h3 - h1;
delta23 = h3 - h2;

% 避免重复计算的中间量
halfBase = 3 * L0 / 2;

%% 构造各段变换矩阵
A1 = eye(4);
A1(3, 4) = h0;

A2 = eye(4);
A2(3, 4) = h1 - h0;

alpha = atan2(delta12, halfBase);
A3 = [1, 0, 0, 0;
      0, cos(alpha), -sin(alpha), -L0 / 2 + L0 * cos(alpha) / 2;
      0, sin(alpha),  cos(alpha),  L0 * sin(alpha) / 2;
      0, 0, 0, 1];

L1 = sqrt(halfBase^2 + delta12^2) - L0 / 2;

%% 摆角相关量
b = cosbeta1hanshu(h1, h2, h3, L0);

triangleEdge = [sqrt(3) * L0 / 2;
                L0 / 2 + L1;
                0];
edgeNorm = norm(triangleEdge);

if edgeNorm < eps
    error('jieA:DegenerateConfiguration', 'Triangle edge norm is near zero.');
end

unitP1 = triangleEdge(1) / edgeNorm;
unitP2 = triangleEdge(2) / edgeNorm;

denominator = sqrt(L0^2 + L1^2 + L0 * L1);
if denominator < eps
    error('jieA:SingularLength', 'Length denominator is near zero.');
end

Lq = (sqrt(3) * L0 * L1 / 2) / denominator;

offsetX = -Lq * (1 - cos(b)) * ((L1 + L0 / 2) / denominator);
offsetY =  Lq * (1 - cos(b)) * ((sqrt(3) * L0 / 2) / denominator);
offsetZ = -Lq * sin(b);

A4 = [unitP1^2 * (1 - cos(b)) + cos(b), unitP1 * unitP2 * (1 - cos(b)), unitP2 * sin(b),  offsetX;
       unitP1 * unitP2 * (1 - cos(b)), unitP2^2 * (1 - cos(b)) + cos(b), -unitP1 * sin(b), offsetY;
      -unitP2 * sin(b),                unitP1 * sin(b),                 cos(b),            offsetZ;
       0,                              0,                              0,                 1];

%% 末端旋转
L13Squared = 3 * L0^2 + delta13^2;
L23Squared = 3 * L0^2 + delta23^2;

[gamma, shiftX, shiftY] = gammafast(h1, h2, h3, L13Squared, L23Squared, L0, L1);

A5 = [cos(gamma), -sin(gamma), 0, shiftX;
      sin(gamma),  cos(gamma), 0, shiftY;
      0,           0,          1, 0;
      0,           0,          0, 1];

%% 组合总变换
T = A1 * A2 * A3 * A4 * A5;
end