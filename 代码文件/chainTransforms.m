function chain = chainTransforms(q, config)
%CHAINTRANSFORMS 计算串联柔性臂各级的齐次变换
%   chain = chainTransforms(q, config) 返回每一级的基座与顶盘在世界坐标
%   下的齐次变换、位置和旋转矩阵。

stageCount = config.stageCount;
numDof = 3 * stageCount;

if numel(q) ~= numDof
	error('chainTransforms:StateDimensionMismatch', ...
		'期望广义坐标维度为 %d，但收到 %d 。', numDof, numel(q));
end

hMatrix = reshape(q, 3, stageCount);

chain(stageCount) = struct('baseTransform', eye(4), ...
	'basePosition', zeros(3, 1), ...
	'baseRotation', eye(3), ...
	'topTransform', eye(4), ...
	'topPosition', zeros(3, 1), ...
	'topRotation', eye(3));

cumulativeTransform = eye(4);

for stageIdx = 1:stageCount
	orientationAngle = config.stage(stageIdx).orientationAngle;
	orientationTransform = rotationAboutZ(orientationAngle);

	baseTransform = cumulativeTransform * orientationTransform;

	h = hMatrix(:, stageIdx);
	geometry = config.stage(stageIdx).geometry;
	localTransform = jieA(h(1), h(2), h(3), geometry.h0, geometry.L0);

	topTransform = baseTransform * localTransform;

	chain(stageIdx).baseTransform = baseTransform;
	chain(stageIdx).basePosition = baseTransform(1:3, 4);
	chain(stageIdx).baseRotation = baseTransform(1:3, 1:3);

	chain(stageIdx).topTransform = topTransform;
	chain(stageIdx).topPosition = topTransform(1:3, 4);
	chain(stageIdx).topRotation = topTransform(1:3, 1:3);

	cumulativeTransform = topTransform;
end
end

function transform = rotationAboutZ(angle)
transform = eye(4);
cosA = cos(angle);
sinA = sin(angle);
transform(1:3, 1:3) = [cosA, -sinA, 0;
					   sinA,  cosA, 0;
					   0,     0,    1];
end
