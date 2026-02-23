function result = simulate_arm(varargin)
%SIMULATE_ARM 运行柔性臂多级动力学仿真，可接收可选参数或结构体
%
%  result = simulate_arm() 使用默认参数运行仿真，并绘制关节与盘心曲线。
%  result = simulate_arm(params) 使用结构体 params 覆盖默认参数字段。
%  result = simulate_arm('stageCount', 6, 'enablePlots', false, ...) 也可
%  采用名称-值对形式传参；可覆盖的字段见 parseSimulationParameters。
%
%  返回 result 结构体，包含时间序列、关节位移/速度以及各级盘心和姿态。
%
%  为支持 Python MATLAB Engine 调用，本函数返回纯数值字段，不依赖句柄。

params = parseSimulationParameters(varargin{:});
params = sanitizeParameters(params);

config = buildSimulationConfig(params);
[timeVector, stateHistory] = ode45(@(t, q) arm_new(t, q, config), ...
	config.timeVector, config.initialState);

stageCount = config.stageCount;
jointDisplacement = stateHistory(:, 1:3 * stageCount);
jointVelocity = stateHistory(:, 3 * stageCount + 1:end);

plateData = quan(stateHistory, config);

if params.enablePlots
	plotJointDisplacement(timeVector, jointDisplacement, config);
	plotPlatePosition(timeVector, plateData, config);
	plotStageTrajectories(timeVector, plateData, config);
end

result = struct( ...
	'time', timeVector, ...
	'jointDisplacement', jointDisplacement, ...
	'jointVelocity', jointVelocity, ...
	'finalPosition', plateData.finalPosition, ...
	'finalVelocity', plateData.finalVelocity, ...
	'orientation', plateData.orientation, ...
	'orientationRate', plateData.orientationRate, ...
	'stagePositions', plateData.stagePositions, ...
	'stageVelocities', plateData.stageVelocities, ...
	'stageCount', stageCount, ...
	'duration', config.duration, ...
	'parameters', params);
end

%% ===================== 参数解析 =====================

function params = parseSimulationParameters(varargin)
	defaults = struct( ...
		'stageCount', 8, ...
		'duration', 10, ...
		'timeStep', 0.01, ...
		'derivativeTimeStep', 1e-4, ...
		'initialHeight', 350e-3, ...
		'initialState', [], ...
		'orientationOffset', pi / 3, ...
		'enablePlots', true, ...
		'customControlProfiles', [], ...
		'customForceMatrix', []); %#ok<*STRNU>

	if nargin == 0
		params = defaults;
		return;
	end

	if nargin == 1 && isstruct(varargin{1})
		overrides = varargin{1};
	else
		if mod(nargin, 2) ~= 0
			error('simulate_arm:InvalidArguments', ...
				'名称-值对需要偶数个输入参数。');
		end
		overrides = struct();
		for idx = 1:2:nargin
			name = varargin{idx};
			value = varargin{idx + 1};
			if ~ischar(name) && ~(isstring(name) && isscalar(name))
				error('simulate_arm:InvalidName', ...
					'名称-值对中的名称必须为字符向量或字符串标量。');
			end
			overrides.(char(name)) = value;
		end
	end

	params = mergeStruct(defaults, overrides);
end

function params = sanitizeParameters(params)
	params.stageCount = max(1, round(params.stageCount));
	params.duration = max(eps, params.duration);
	params.timeStep = max(eps, params.timeStep);
	params.derivativeTimeStep = max(eps, params.derivativeTimeStep);
	params.initialHeight = double(params.initialHeight);

	if isempty(params.initialState)
		params.initialState = [];
	else
		params.initialState = double(params.initialState(:));
	end

	params.orientationOffset = double(params.orientationOffset);
	params.enablePlots = logical(params.enablePlots);

	if isempty(params.customForceMatrix)
		params.customForceMatrix = [];
	elseif isnumeric(params.customForceMatrix)
		params.customForceMatrix = double(params.customForceMatrix);
		if size(params.customForceMatrix, 2) ~= 3
			error('simulate_arm:InvalidCustomForceMatrix', ...
				'customForceMatrix 必须为 [stageCount × 3] 的矩阵。');
		end
	elseif iscell(params.customForceMatrix)
		params.customForceMatrix = reshape(params.customForceMatrix, [], 1);
	elseif isstring(params.customForceMatrix)
		params.customForceMatrix = cellstr(params.customForceMatrix(:));
	elseif ischar(params.customForceMatrix)
		params.customForceMatrix = {params.customForceMatrix};
	else
		error('simulate_arm:InvalidCustomForceMatrix', ...
			'customForceMatrix 必须为数值矩阵或包含表达式的单元数组。');
	end

	if isempty(params.customControlProfiles)
		params.customControlProfiles = [];
	elseif ~iscell(params.customControlProfiles)
		error('simulate_arm:InvalidControlProfiles', ...
			'customControlProfiles 必须为 Cell 数组。');
	end
end

function merged = mergeStruct(defaults, overrides)
	merged = defaults;
	if isempty(overrides)
		return;
	end

	overrideFields = fieldnames(overrides);
	for idx = 1:numel(overrideFields)
		field = overrideFields{idx};
		merged.(field) = overrides.(field);
	end
end

%% ===================== 配置与仿真 =====================

function config = buildSimulationConfig(params)
	stageCount = params.stageCount;

	config.stageCount = stageCount;
	config.orientationOffset = params.orientationOffset;
	config.timeStep = params.timeStep;
	config.duration = params.duration;
	config.timeVector = 0:config.timeStep:config.duration;
	config.derivativeTimeStep = params.derivativeTimeStep;

	stageTemplate = createStageTemplate(params.initialHeight);
	config.stage = repmat(stageTemplate, stageCount, 1);

	defaultInitialH = params.initialHeight;
	defaultState = [repmat(defaultInitialH, 3 * stageCount, 1);
				 zeros(3 * stageCount, 1)];

	if ~isempty(params.initialState)
		if numel(params.initialState) ~= 6 * stageCount
			error('simulate_arm:InitialStateDimensionMismatch', ...
				'initialState 维度应为 %d，但收到 %d。', 6 * stageCount, numel(params.initialState));
		end
		config.initialState = params.initialState;
	else
		config.initialState = defaultState;
	end

	controlProfiles = createStageControlProfiles(stageCount, params.customControlProfiles, params.customForceMatrix);

	for stage = 1:stageCount
		config.stage(stage).geometry.h0 = defaultInitialH;
		config.stage(stage).orientationAngle = (stage - 1) * config.orientationOffset;
		config.stage(stage).control = mergeControlTemplate( ...
			config.stage(stage).control, controlProfiles{stage});
	end
end

function template = createStageTemplate(initialHeight)
	template.geometry = struct( ...
		'h0', initialHeight, ...
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

function controlProfiles = createStageControlProfiles(stageCount, customProfiles, customForceMatrix)
	baseTarget = [-0.08; 0.1; 5.0];
	baseAmplitudes = [0.12; 0.07; 0.05];
	baseFrequencies = [pi; 0.7 * pi; 0.45 * pi];
	basePhases = [0; pi / 6; pi / 4];

	controlProfiles = cell(stageCount, 1);

	for stage = 1:stageCount
		profile = struct('target', [], 'customInput', [], 'externalDirection', []);
		constantForce = [];
		expressionInput = [];
		if ~isempty(customForceMatrix)
			if isnumeric(customForceMatrix) && size(customForceMatrix, 1) >= stage
				row = customForceMatrix(stage, :);
				if any(isfinite(row))
					constantForce = row(:);
				end
			elseif iscell(customForceMatrix) && numel(customForceMatrix) >= stage
				entry = customForceMatrix{stage};
				[constantForce, expressionInput] = interpretCustomForceEntry(entry);
			end
		end

		if ~isempty(customProfiles) && numel(customProfiles) >= stage && ~isempty(customProfiles{stage})
			customProfile = customProfiles{stage};
			overrideFields = fieldnames(customProfile);
			for k = 1:numel(overrideFields)
				field = overrideFields{k};
				profile.(field) = customProfile.(field);
			end
		end

		if isempty(profile.target)
			targetOffset = [0.005 * (stage - 1);
						 0.004 * (stage - 1);
						 0.06 * (stage - 1)];
			profile.target = baseTarget + targetOffset;
		end

		if isempty(profile.externalDirection)
			profile.externalDirection = [1; 0; 0];
		end

		if isempty(profile.customInput)
			amplitudeVec = baseAmplitudes .* (1 + 0.08 * (stage - 1));
			frequencyVec = baseFrequencies .* (1 + 0.05 * (stage - 1));
			phaseVec = basePhases + (stage - 1) * (pi / 10);
			amplitudeCopy = amplitudeVec;
			frequencyCopy = frequencyVec;
			phaseCopy = phaseVec;
			profile.customInput = @(time, stageState) ...
				amplitudeCopy .* sin(frequencyCopy * time + phaseCopy);
		end

		if ~isempty(expressionInput)
			profile.customInput = expressionInput;
		elseif ~isempty(constantForce)
			profile.customInput = constantForce;
		end

		profile.customInput = finalizeCustomInput(profile.customInput);

		controlProfiles{stage} = profile;
	end
end

function control = mergeControlTemplate(baseControl, override)
	control = baseControl;
	overrideFields = fieldnames(override);
	for idx = 1:numel(overrideFields)
		field = overrideFields{idx};
		if ~isempty(override.(field))
			control.(field) = override.(field);
		end
	end
end

%% ===================== 绘图辅助 =====================

function plotJointDisplacement(timeVector, jointDisplacement, config)
	stageCount = config.stageCount;
	positionMm = jointDisplacement * 1e3;

	figure('Name', '柔性臂动力学仿真结果 - 关节', 'Color', 'w');
	tiledlayout(stageCount, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
	for stage = 1:stageCount
		for rod = 1:3
			nexttile;
			columnIndex = (stage - 1) * 3 + rod;
			plot(timeVector, positionMm(:, columnIndex), 'LineWidth', 2);
			title(sprintf('第%d级第%d根杆位置 (mm)', stage, rod));
			xlabel('时间 (s)');
			ylabel('位移 (mm)');
			grid on;
		end
	end
end

function plotPlatePosition(timeVector, plateData, config)
	stageCount = config.stageCount;
	positionMm = plateData.stagePositions * 1e3;
	colorOrder = lines(stageCount);

	figure('Name', '柔性臂动力学仿真结果 - 盘心', 'Color', 'w');
	tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
	componentLabels = {'x', 'y', 'z'};
	for axisIdx = 1:3
		nexttile;
		hold on;
		for stage = 1:stageCount
			series = squeeze(positionMm(:, axisIdx, stage));
			lineWidth = 1.6 + 0.6 * (stage == stageCount);
			plot(timeVector, series, 'LineWidth', lineWidth, 'Color', colorOrder(stage, :));
		end
		hold off;
		xlabel('时间 (s)');
		ylabel('位移 (mm)');
		title(sprintf('%s 方向盘心位置变化', componentLabels{axisIdx}));
		legend(arrayfun(@(s) sprintf('第%d级盘心', s), 1:stageCount, 'UniformOutput', false), ...
			'Location', 'best');
		grid on;
	end
end

function plotStageTrajectories(timeVector, plateData, config)
	stageCount = config.stageCount;
	finalPosition = plateData.finalPosition;

	figure('Name', '盘心三维轨迹', 'Color', 'w');
	ax = axes('Parent', gcf);
	hold(ax, 'on');
	colorOrder = lines(stageCount);
	for stage = 1:stageCount
		trajectory = squeeze(plateData.stagePositions(:, :, stage));
		plot3(ax, trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), ...
			'LineWidth', 1.6, 'Color', colorOrder(stage, :));
	end
	scatter3(ax, finalPosition(:, 1), finalPosition(:, 2), finalPosition(:, 3), 18, timeVector, ...
		'filled', 'DisplayName', sprintf('第%d级散点', stageCount));
	xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)'); zlabel(ax, 'Z (m)');
	title(ax, '各级盘心位置轨迹'); grid(ax, 'on'); axis(ax, 'equal');
	legendEntries = arrayfun(@(s) sprintf('第%d级轨迹', s), 1:stageCount, 'UniformOutput', false);
	legend(ax, [legendEntries, {sprintf('第%d级散点', stageCount)}], 'Location', 'best');
	hold(ax, 'off');
end

function [constantForce, expressionInput] = interpretCustomForceEntry(entry)
	constantForce = [];
	expressionInput = [];

	if isempty(entry)
		return;
	end

	if iscell(entry)
		if numel(entry) == 3
			allNumeric = all(cellfun(@isnumeric, entry));
			allString = all(cellfun(@(v) ischar(v) || isstring(v), entry));
			if allNumeric
				constantForce = cellfun(@double, entry(:));
			elseif allString
				expressionInput = composeVectorExpression(entry(:));
			end
		end
	elseif isnumeric(entry)
		vector = double(entry(:));
		if numel(vector) == 3
			constantForce = vector;
		end
	elseif isstring(entry)
		if any(strlength(entry) > 0)
			expressionInput = char(entry(1));
		end
	elseif ischar(entry)
		trimmed = strtrim(entry);
		if ~isempty(trimmed)
			expressionInput = trimmed;
		end
	end
end

function expression = composeVectorExpression(components)
	stringComponents = cellfun(@(c) strtrim(char(c)), components, ...
		'UniformOutput', false);
	expression = ['[', strjoin(stringComponents, '; '), ']'];
end

function handle = finalizeCustomInput(customInput)
	if isempty(customInput)
		handle = [];
		return;
	end

	if isa(customInput, 'function_handle')
		handle = customInput;
		return;
	end

	if isnumeric(customInput)
		constantVector = validateGeneralizedVector(customInput, 'numeric customInput');
		handle = @(time, stageState) constantVector;
		return;
	end

	if iscell(customInput)
		numericCandidate = cellfun(@(v) double(v), customInput, 'UniformOutput', false);
		numericCandidate = vertcat(numericCandidate{:});
		handle = finalizeCustomInput(numericCandidate);
		return;
	end

	if ischar(customInput) || isstring(customInput)
		handle = createExpressionHandle(customInput);
		return;
	end

	error('simulate_arm:UnsupportedCustomInput', '无法识别 customInput 的类型。');
end

function handle = createExpressionHandle(expression)
	expressionChar = char(expression);
	handle = @(time, stageState) evaluateCustomExpression(expressionChar, time, stageState);
end

function vector = evaluateCustomExpression(expression, time, stageState)
	
	t = time; %#ok<NASGU>
	h = stageState.h; %#ok<NASGU>
	hd = stageState.hd; %#ok<NASGU>
	vector = eval(expression);
	vector = validateGeneralizedVector(vector, 'expression customInput');
end

function vector = validateGeneralizedVector(value, contextLabel)
	vector = double(value(:));
	if numel(vector) ~= 3 || any(~isfinite(vector))
		error('simulate_arm:InvalidCustomInputVector', ...
			'%s 必须生成 3×1 的有限值向量。', contextLabel);
	end
end
