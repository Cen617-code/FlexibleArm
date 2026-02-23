
function dq = arm_new(t, qh, config)
%ARM_NEW 多级柔性臂动力学模型（基于拉格朗日法）
%   dq = arm_new(t, qh, config) 返回状态向量 qh 的导数，用于 ODE 积分。
%   qh = [h; h_dot]，其中前 3*stageCount 个分量为各级广义坐标，后 3*stageCount
%   个分量为其一阶导。

if nargin < 3 || isempty(config)
	config = buildDefaultArmConfig();
end

stageCount = config.stageCount;
stateDimension = 6 * stageCount;
if numel(qh) ~= stateDimension
	error('arm_new:StateDimensionMismatch', ...
		'期望状态维度为 %d，但收到 %d 。', stateDimension, numel(qh));
end

%% 状态拆分
stageDisplacement = reshape(qh(1:3 * stageCount), 3, stageCount);
stageVelocity = reshape(qh(3 * stageCount + 1:end), 3, stageCount);

%% 数值配置
fd = loadFiniteDifferenceConfig(config);

%% 装配动力学系统
system = assembleLagrangeSystem(t, stageDisplacement, stageVelocity, config, fd);

massMatrix = system.massMatrix;
generalizedForce = system.generalizedForce;
constraintJacobian = system.constraintJacobian;
constraintDrift = system.constraintDrift;

numDof = size(massMatrix, 1);

massMatrix = regularizeSymmetricMatrix(massMatrix, 1e-8);
generalizedForce = sanitizeVector(generalizedForce);
constraintJacobian = sanitizeMatrix(constraintJacobian);
constraintDrift = sanitizeVector(constraintDrift);

if isempty(constraintJacobian)
	accelerationVector = solveLinearSystem(massMatrix, generalizedForce);
else
	zeroBlock = zeros(size(constraintJacobian, 1));
	augmentedMatrix = [massMatrix, constraintJacobian';
					  constraintJacobian, zeroBlock];
	augmentedMatrix = regularizeSymmetricMatrix(augmentedMatrix, 1e-10);
	rhs = [generalizedForce;
		  -constraintDrift];
	solution = solveLinearSystem(augmentedMatrix, rhs);
	accelerationVector = solution(1:numDof);
end

stageAcceleration = reshape(accelerationVector, 3, stageCount);

%% 组装导数
dq = [stageVelocity(:);
	  stageAcceleration(:)];
end

%% ===================== 辅助函数 =====================

function system = assembleLagrangeSystem(t, stageDisplacement, stageVelocity, config, fd)
stageCount = config.stageCount;
numDof = 3 * stageCount;

q = stageDisplacement(:);
qdot = stageVelocity(:);

chainNominal = chainTransforms(q, config);

basePositions = zeros(3, stageCount);
topPositions = zeros(3, stageCount);
topRotations = zeros(3, 3, stageCount);

for stageIdx = 1:stageCount
	basePositions(:, stageIdx) = chainNominal(stageIdx).basePosition;
	topPositions(:, stageIdx) = chainNominal(stageIdx).topPosition;
	topRotations(:, :, stageIdx) = chainNominal(stageIdx).topRotation;
end

delta = fd.spatialStep;

baseJacobians = zeros(3, numDof, stageCount);
topJacobians = zeros(3, numDof, stageCount);
rotationJacobians = zeros(3, numDof, stageCount);

for dofIdx = 1:numDof
	perturbedQ = q;
	perturbedQ(dofIdx) = perturbedQ(dofIdx) + delta;

	chainPerturbed = chainTransforms(perturbedQ, config);

	for stageIdx = 1:stageCount
		baseJacobians(:, dofIdx, stageIdx) = ...
			(chainPerturbed(stageIdx).basePosition - basePositions(:, stageIdx)) / delta;

		topJacobians(:, dofIdx, stageIdx) = ...
			(chainPerturbed(stageIdx).topPosition - topPositions(:, stageIdx)) / delta;

		rotationDifference = topRotations(:, :, stageIdx)' * chainPerturbed(stageIdx).topRotation;
		rotationJacobians(:, dofIdx, stageIdx) = logSO3(rotationDifference) / delta;
	end
end

massMatrix = zeros(numDof);
generalizedForce = zeros(numDof, 1);

for stageIdx = 1:stageCount
	stageRange = (stageIdx - 1) * 3 + (1:3);
	massParams = config.stage(stageIdx).mass;
	controlParams = config.stage(stageIdx).control;

	Jb = squeeze(baseJacobians(:, :, stageIdx));
	Jt = squeeze(topJacobians(:, :, stageIdx));
	Jr = squeeze(rotationJacobians(:, :, stageIdx));

	Jb = sanitizeMatrix(Jb);
	Jt = sanitizeMatrix(Jt);
	Jr = sanitizeMatrix(Jr);

	massMatrix = massMatrix + massParams.body * (Jb' * Jb) + massParams.platform * (Jt' * Jt);

	rotationWorld = topRotations(:, :, stageIdx);
	inertiaWorld = rotationWorld * massParams.inertia * rotationWorld';
	massMatrix = massMatrix + Jr' * inertiaWorld * Jr;

	position = topPositions(:, stageIdx);
	positionJacobianLocal = sanitizeMatrix(Jt(:, stageRange));

	stageState.h = stageDisplacement(:, stageIdx);
	stageState.hd = stageVelocity(:, stageIdx);

	controlForce = computeControlForce(t, stageState, position, positionJacobianLocal, controlParams);

	generalizedForce(stageRange) = generalizedForce(stageRange) + controlForce;
end

% 构造逐级连接的约束雅可比
if stageCount > 1
	constraintJacobian = zeros(3 * (stageCount - 1), numDof);
	constraintDrift = zeros(3 * (stageCount - 1), 1);

	qPredicted = q + qdot * fd.timeStep;
	chainPredicted = chainTransforms(qPredicted, config);

	for stageIdx = 1:stageCount - 1
		rowRange = (stageIdx - 1) * 3 + (1:3);

		JTopCurrent = squeeze(topJacobians(:, :, stageIdx));
		JBaseNext = squeeze(baseJacobians(:, :, stageIdx + 1));
		constraintJacobian(rowRange, :) = sanitizeMatrix(JTopCurrent - JBaseNext);

		phiCurrent = topPositions(:, stageIdx) - basePositions(:, stageIdx + 1);
		phiPredicted = chainPredicted(stageIdx).topPosition - chainPredicted(stageIdx + 1).basePosition;
		constraintDrift(rowRange) = sanitizeVector((phiPredicted - phiCurrent) / fd.timeStep);
	end
else
	constraintJacobian = [];
	constraintDrift = [];
end

% 轻微对称化以消除数值误差
massMatrix = regularizeSymmetricMatrix(massMatrix, 1e-12);

system.massMatrix = massMatrix;
system.generalizedForce = generalizedForce;
system.constraintJacobian = constraintJacobian;
system.constraintDrift = constraintDrift;
end

function defaultConfig = buildDefaultArmConfig()
	defaultConfig.stageCount = 1;
	defaultConfig.derivativeTimeStep = 1e-6;
	stageTemplate = createStageTemplate();
	defaultConfig.stage = stageTemplate;
end

function template = createStageTemplate()
	template.geometry = struct( ...
		'h0', 350e-3, ...
		'L0', 55e-3);
	template.mass = struct( ...
		'body', 1.774, ...
		'platform', 9.08, ...
		'inertia', diag([36825.653e-6, 36825.653e-6, 72966.445e-6]));
	template.control = struct( ...
		'gain', 0, ...
		'stiffness', -12.2, ...
		'damping', 18.3, ...
		'target', [-0.08; 0.1; 5.0], ...
		'externalAmplitude', 0.1045, ...
		'externalFrequency', pi, ...
		'externalPhase', 0, ...
		'externalDirection', [1; 0; 0], ...
		'customInput', []);
	template.external = struct('distributedForce', [0; 0; 0]);
	template.orientationAngle = 0;
end

function fd = loadFiniteDifferenceConfig(config)
	fd.spatialStep = 1e-6;
	if isfield(config, 'derivativeTimeStep') && ~isempty(config.derivativeTimeStep)
		fd.timeStep = config.derivativeTimeStep;
	else
		fd.timeStep = 1e-6;
	end
end

function controlForce = computeControlForce(t, state, position, positionJacobian, control)
	displacement = position - control.target;

	positionJacobian = sanitizeMatrix(positionJacobian);
	rcondValue = rcond(positionJacobian);

	if ~isfinite(rcondValue) || rcondValue < 1e-8
		generalizedDisplacement = pinv(positionJacobian, 1e-8) * displacement;
	else
		generalizedDisplacement = positionJacobian \ displacement;
	end

	baseForce = control.gain * (control.stiffness * generalizedDisplacement - control.damping * state.hd);

	if isfield(control, 'customInput') && ~isempty(control.customInput)
		externalGeneralized = control.customInput(t, state);
		externalGeneralized = externalGeneralized(:);
		validateExternalInput(externalGeneralized, 'customInput');
	elseif isfield(control, 'externalDirection') && ~isempty(control.externalDirection)
		direction = control.externalDirection(:);
		if numel(direction) ~= 3
			error('arm_new:InvalidExternalDirection', 'control.externalDirection 必须为 3×1 向量。');
		end
		scalarInput = control.externalAmplitude * sin(control.externalFrequency * t + control.externalPhase);
		externalGeneralized = direction * scalarInput;
	else
		scalarInput = control.externalAmplitude * sin(control.externalFrequency * t + control.externalPhase);
		externalGeneralized = [scalarInput; 0; 0];
	end

	controlForce = baseForce + externalGeneralized;
end

function validateExternalInput(vector, label)
	if numel(vector) ~= 3 || ~isvector(vector)
		error('arm_new:ExternalInputDimensionError', ...
			'%s 必须返回 3×1 向量。', label);
	end
end

function phi = logSO3(rotationMatrix)
	cosTheta = (trace(rotationMatrix) - 1) / 2;
	cosTheta = max(min(cosTheta, 1), -1);
	theta = acos(cosTheta);

	if theta < 1e-8
		phi = 0.5 * [rotationMatrix(3, 2) - rotationMatrix(2, 3);
				 rotationMatrix(1, 3) - rotationMatrix(3, 1);
				 rotationMatrix(2, 1) - rotationMatrix(1, 2)];
	else
		phi = theta / (2 * sin(theta)) * [rotationMatrix(3, 2) - rotationMatrix(2, 3);
			 rotationMatrix(1, 3) - rotationMatrix(3, 1);
			 rotationMatrix(2, 1) - rotationMatrix(1, 2)];
	end
end

function matrix = sanitizeMatrix(matrix)
    invalidMask = ~isfinite(matrix);
    if any(invalidMask, 'all')
        matrix(invalidMask) = 0;
    end
end

function vector = sanitizeVector(vector)
	if isempty(vector)
		return;
	end
	invalidMask = ~isfinite(vector);
	if any(invalidMask)
		vector(invalidMask) = 0;
	end
end

function matrix = regularizeSymmetricMatrix(matrix, epsilon)
	if isempty(matrix)
		return;
	end
	matrix = (matrix + matrix') / 2;
	if epsilon > 0
		matrix = matrix + eye(size(matrix, 1)) * epsilon;
	end
end

function solution = solveLinearSystem(matrix, rhs)
	rcondValue = rcond(matrix);
	if ~isfinite(rcondValue) || rcondValue < 1e-10
		solution = pinv(matrix, 1e-10) * rhs;
	else
		solution = matrix \ rhs;
	end
end