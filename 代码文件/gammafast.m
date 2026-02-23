function [gamma, offsetX, offsetY] = gammafast(h1, h2, h3, ~, ~, L0, L1)
%GAMMAFAST 计算柔性臂末端转角及偏移
%   [gamma, offsetX, offsetY] = gammafast(...) 返回末端平面内的旋转角和
%   与参考点的平面偏移，输入为三广义位移及导出的几何量。

%% 几何前置
if any([L0, L1] <= 0)
    error('gammafast:InvalidGeometry', 'L0 and L1 must be positive.');
end

[x13, y13] = computeP1P3Projection(h1, h2, h3, L0);

%% 通过牛顿迭代求平面偏移
[offsetX, offsetY] = niu([0, 0], L0, L1, x13, y13);

%% 计算末端旋转角
vectorOcP2 = [0; L1; 0];
vectorOcpP2 = [-offsetX; L1 - offsetY; 0];

crossVector = cross(vectorOcP2, vectorOcpP2);
dotValue = dot(vectorOcP2, vectorOcpP2);

normProduct = norm(vectorOcP2) * norm(vectorOcpP2);
if normProduct < eps
    error('gammafast:DegenerateConfiguration', 'Vectors for gamma computation are nearly zero.');
end

cosGamma = max(min(dotValue / normProduct, 1), -1);

if crossVector(3) >= 0
    gamma = acos(cosGamma);
else
    gamma = -acos(cosGamma);
end

gamma = real(gamma);
end

%% ===================== 辅助函数 =====================

function [x13, y13] = computeP1P3Projection(h1, h2, h3, L0)
    sqrt3 = sqrt(3);

    term1 = 2 * h1^2 - 2 * h1 * h2 - 2 * h1 * h3 + 2 * h2 * h3;
    term2 = 3 * L0^2;

    delta12 = h1 - h2;

    sqrtAArg = 4 * delta12^2 + 9 * L0^2;
    sqrtA = sqrt(max(sqrtAArg, 0));

    sqrtBArg = L0^2 * (4 * h1^2 - 4 * h1 * h2 - 4 * h1 * h3 + 4 * h2^2 - 4 * h2 * h3 + 4 * h3^2 + 9 * L0^2);
    sqrtB = sqrt(max(sqrtBArg, 0));

    denominator = 4 * ((h1 - h2)^2 + 3 * L0^2) * L0;
    if abs(denominator) < eps
        error('gammafast:SingularProjection', 'Denominator in projection computation is near zero.');
    end

    x13 = sqrt3 * (term1 * L0^2 + term2 * L0^2 + sqrtA * sqrtB) / denominator;

    numeratorY = 2 * sqrtA * h1^2 - 2 * sqrtA * h1 * h2 - 2 * sqrtA * h1 * h3 + 2 * sqrtA * h2 * h3 + 3 * sqrtA * L0^2 - 3 * sqrtB;
    denominatorY = 4 * ((h1 - h2)^2 + 3 * L0^2);

    if abs(denominatorY) < eps
        error('gammafast:SingularProjection', 'Denominator in Y projection is near zero.');
    end

    y13 = numeratorY / denominatorY;
end
% syms a;
% OcOcp = [x;y];%所有都在Oc系下
% OOc = [sqrt(3)*L0/3+sqrt(3)*L1/6;-L1/2];
% OOcz = [L1/2;sqrt(3)*L0/3+sqrt(3)*L1/6];
% OOcp = OOc+OcOcp;
% OOc=OOc./sqrt(sum(OOc.^2));
% OOcz=OOcz./sqrt(sum(OOcz.^2));
% OOcp=double(OOcp./sqrt(sum(OOcp.^2)));
% f1 = cos(a)*OOc(1)+sin(a)*OOcz(1)-OOcp(1);
% f2 = cos(a)*OOc(2)+sin(a)*OOcz(2)-OOcp(2);
% r1 =double(solve(f1,"real",true));
% r2=double(solve(f2,"real",true));
% %寻找公共解
% rp1=[r1;r2];
% Rp1=sort(rp1);
% dr=diff(Rp1);
% dr=abs(dr);
% [~,ii]=sort(dr);
% r=(Rp1(ii(1))+Rp1(ii(1)+1))/2;
% %最后这个为啥就无解呢啊啊啊啊啊啊啊啊
%把x和y分别带进去求解算的解完全不一样啊
% %%图像
% aa=0:0.01:8;
% yy1=cos(aa).*OOc(1)+sin(aa).*OOcz(1)-OOcp(1);
% yy2=cos(aa).*OOc(2)+sin(aa).*OOcz(2)-OOcp(2);
% plot(aa,yy1)
% hold on
% plot(aa,yy2)
% plot(aa,zeros(length(aa),1))
% hold off