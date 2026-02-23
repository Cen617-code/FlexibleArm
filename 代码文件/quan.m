function varargout = quan(stateHistory, config)
%QUAN 计算柔性臂仿真数据中的多级盘心与姿态派生量
%   当仅请求一个输出时，返回结构体，包含各级盘心位置/速度以及最终级的
%   姿态信息。保留旧版 12 个输出参数的兼容接口。

if nargin < 2 || isempty(config)
    config = buildDefaultQuanConfig();
end

stageCount = config.stageCount;
expectedColumns = stageCount * 6;
if size(stateHistory, 2) ~= expectedColumns
    error('quan:StateDimensionMismatch', ...
        '期望输入包含 %d 列 (3h + 3hd)，但收到 %d 列。', expectedColumns, size(stateHistory, 2));
end

if isfield(config, 'derivativeTimeStep') && ~isempty(config.derivativeTimeStep)
    deltaT = config.derivativeTimeStep;
elseif isfield(config, 'timeStep') && ~isempty(config.timeStep)
    deltaT = min(config.timeStep, 1e-4);
else
    deltaT = 1e-4;
end

if deltaT <= 0
    error('quan:InvalidTimeStep', '有限差分时间步长必须为正。');
end

numSamples = size(stateHistory, 1);

finalPosition = zeros(numSamples, 3);
finalVelocity = zeros(numSamples, 3);
orientationVector = zeros(numSamples, 3);
orientationRate = zeros(numSamples, 3);
stagePositions = zeros(numSamples, 3, stageCount);
stageVelocities = zeros(numSamples, 3, stageCount);

for idx = 1:numSamples
    stageH = reshape(stateHistory(idx, 1:3 * stageCount), 3, stageCount);
    stageHd = reshape(stateHistory(idx, 3 * stageCount + 1:end), 3, stageCount);

    predictedH = stageH + stageHd * deltaT;

    [currentPositions, currentRotations] = computeStageChain(stageH, config);
    [predPositions, ~] = computeStageChain(predictedH, config);

    for stageIdx = 1:stageCount
        positionCurrent = currentPositions(:, stageIdx);
        positionNext = predPositions(:, stageIdx);

        stagePositions(idx, :, stageIdx) = positionCurrent.';
        stageVelocities(idx, :, stageIdx) = ((positionNext - positionCurrent) / deltaT).';
    end

    finalCurrent = currentPositions(:, stageCount);
    finalNext = predPositions(:, stageCount);

    finalPosition(idx, :) = finalCurrent.';
    finalVelocity(idx, :) = ((finalNext - finalCurrent) / deltaT).';

    finalRotation = currentRotations(:, :, stageCount);
    rotationWorldToBody = finalRotation.';

    geometry = config.stage(stageCount).geometry;

    axisA = normalizeVector(rotationWorldToBody * [1; 0; 0]);
    axisB = normalizeVector(rotationWorldToBody * [sqrt(3) * geometry.L0 / 2;
                              3 * geometry.L0 / 2;
                              stageH(2, stageCount) - stageH(1, stageCount)]);
    axisR = [0; 0; 1];
    axisTransform = [axisA, axisB, axisR];

    anglesCurrent = evaluateAngles(stageH(:, stageCount), geometry);
    anglesNext = evaluateAngles(predictedH(:, stageCount), geometry);

    orientationNow = axisTransform * anglesCurrent;
    orientationNext = axisTransform * anglesNext;

    orientationVector(idx, :) = orientationNow.';
    orientationRate(idx, :) = ((orientationNext - orientationNow) / deltaT).';
end

dataStruct = struct( ...
    'finalPosition', finalPosition, ...
    'finalVelocity', finalVelocity, ...
    'orientation', orientationVector, ...
    'orientationRate', orientationRate, ...
    'stagePositions', stagePositions, ...
    'stageVelocities', stageVelocities, ...
    'stageCount', stageCount);

if nargout <= 1
    varargout{1} = dataStruct;
    return;
end

varargout{1} = finalPosition(:, 1);
varargout{2} = finalPosition(:, 2);
varargout{3} = finalPosition(:, 3);
varargout{4} = orientationVector(:, 1);
varargout{5} = orientationVector(:, 2);
varargout{6} = orientationVector(:, 3);
varargout{7} = finalVelocity(:, 1);
varargout{8} = finalVelocity(:, 2);
varargout{9} = finalVelocity(:, 3);
varargout{10} = orientationRate(:, 1);
varargout{11} = orientationRate(:, 2);
varargout{12} = orientationRate(:, 3);
if nargout >= 13
    varargout{13} = dataStruct;
end
end

%% ===================== 辅助函数 =====================

function [positions, rotations] = computeStageChain(hStages, config)
    stageCount = config.stageCount;
    positions = zeros(3, stageCount);
    rotations = zeros(3, 3, stageCount);
    transformAccum = eye(4);

    for stageIdx = 1:stageCount
        geometry = config.stage(stageIdx).geometry;
        orientationAngle = 0;
        if isfield(config.stage(stageIdx), 'orientationAngle') && ~isempty(config.stage(stageIdx).orientationAngle)
            orientationAngle = config.stage(stageIdx).orientationAngle;
        end

        localTransform = jieA(hStages(1, stageIdx), hStages(2, stageIdx), hStages(3, stageIdx), geometry.h0, geometry.L0);
        rotationTransform = makeRotationAboutZ(orientationAngle);
        stageTransform = rotationTransform * localTransform;

        transformAccum = transformAccum * stageTransform;
        positions(:, stageIdx) = transformAccum(1:3, 4);
        rotations(:, :, stageIdx) = transformAccum(1:3, 1:3);
    end
end

function transform = makeRotationAboutZ(angle)
    c = cos(angle);
    s = sin(angle);
    transform = eye(4);
    transform(1:3, 1:3) = [c, -s, 0;
                             s,  c, 0;
                             0,  0, 1];
end

function angles = evaluateAngles(h, geometry)
    [a, b, r] = jieabr(h(1), h(2), h(3), geometry.L0);
    angles = [a; b; r];
end

function vector = normalizeVector(vector)
    magnitude = norm(vector);
    if magnitude < eps
        return;
    end
    vector = vector / magnitude;
end

function defaultConfig = buildDefaultQuanConfig()
    defaultConfig.stageCount = 1;
    defaultConfig.derivativeTimeStep = 1e-4;
    stageTemplate = createStageTemplate();
    stageTemplate.geometry.h0 = 350e-3;
    defaultConfig.stage = stageTemplate;
end

function template = createStageTemplate()
    template.geometry = struct( ...
        'h0', 0, ...
        'L0', 55e-3);
    template.mass = struct();
    template.control = struct();
    template.external = struct();
    template.orientationAngle = 0;
end