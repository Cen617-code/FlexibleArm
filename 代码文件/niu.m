function [x, y] = niu(initialGuess, L0, L1, x13, y13)
%NIU 通过牛顿迭代求解柔性臂盘心在 Oc 坐标系下的平面坐标

%% 常量与初始化
tolerance = 1e-20;
finiteStep = 1e-5;
maxIterations = 200;

iterationState = initialGuess(:);

vectorP1P3 = [x13; y13];
vectorP1P2 = [sqrt(3) * L0 / 2;
              L0 / 2 + L1];

%% 迭代求解
for iteration = 1:maxIterations
    [residual, jacobian] = evaluateResidualAndJacobian(iterationState, vectorP1P2, vectorP1P3, L0, finiteStep);

    if norm(residual)^2 <= tolerance
        break;
    end

    rcondValue = rcond(jacobian);
    if ~isfinite(rcondValue) || rcondValue < 1e-10
        lambda = 1e-8;
        regularizedJacobian = jacobian + lambda * eye(2);
        delta = -pinv(regularizedJacobian) * residual;
    else
        delta = -jacobian \ residual;
    end
    iterationState = iterationState + delta;
end

if norm(residual)^2 > tolerance
    warning('niu:MaxIterationsReached', 'Newton solver reached max iterations without convergence.');
end

x = iterationState(1);
y = iterationState(2);
end

%% ===================== 辅助函数 =====================

function [residualVector, jacobianMatrix] = evaluateResidualAndJacobian(point, vectorP1P2, vectorP1P3, L0, finiteStep)
    baseResidual = evaluateResidual(point, vectorP1P2, vectorP1P3, L0);

    residualVector = baseResidual;
    jacobianMatrix = zeros(2, 2);

    for idx = 1:2
        perturbation = zeros(2, 1);
        perturbation(idx) = finiteStep;

        perturbedPoint = point + perturbation;
        perturbedResidual = evaluateResidual(perturbedPoint, vectorP1P2, vectorP1P3, L0);

        jacobianMatrix(:, idx) = (perturbedResidual - baseResidual) / finiteStep;
    end
end

function residualVector = evaluateResidual(point, vectorP1P2, vectorP1P3, L0)
    ocpP1 = [-sqrt(3) * L0 / 2 - point(1);
             -L0 / 2 - point(2)];

    termA = ocpP1 + vectorP1P3;
    termB = ocpP1 + vectorP1P2;

    residual1 = 0.5 * norm(termA) * norm(termB) + dot(termA, termB);
    residual2 = 0.5 * norm(ocpP1) * norm(termB) + dot(ocpP1, termB);

    residualVector = [residual1; residual2];
end