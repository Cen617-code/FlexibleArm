# Python 图形界面：柔性臂仿真控制台

该目录提供一个基于 `tkinter` 的轻量 GUI，用于设置柔性臂动力学仿真的关键参数，并调用 MATLAB 脚本运行 `simulate_arm.m`。

## 依赖与准备

1. **Python 3.8+**（Windows 自带 `tkinter`，无需额外安装）。
2. **MATLAB Engine for Python**：
   - 在 MATLAB 安装目录下执行：
     ```powershell
     cd "<MATLABROOT>\extern\engines\python"
     python setup.py install
     ```
   - 安装成功后，Python 中可 `import matlab.engine`。
3. 将代码仓库解压或克隆至本地，保持 MATLAB `.m` 文件与 `python_ui` 目录在同一层级。

> 若尚未配置 MATLAB 命令行，可直接依赖 MATLAB Engine，无需将 `matlab.exe` 加入 PATH。

## 使用方法

```powershell
cd path\to\柔性臂动力学仿真系列函数\二层\python_ui
python app.py
```

界面提供以下可调参数：

| 参数 | 说明 |
| ---- | ---- |
| Stage Count | 串联级数，默认 8 |
| 仿真时长 | `simulate_arm` 中的 `duration`，单位秒 |
| 积分步长 dt | ODE45 时间步配置，默认 0.01 s |
| 有限差分步长 | 传递给 `arm_new` 的导数步长，默认 1e-4 s |
| 初始高度 h0 | 各级 `geometry.h0` 默认值，单位米 |
| 显示 MATLAB 绘图 | 勾选后沿用 MATLAB Figure 输出；未勾选时仅计算数据 |
| 各杆力表达式 f_z(t) | 为每级三根杆指定沿局部 Z 轴的力函数，可使用变量 `t` 及 MATLAB 表达式；留空表示沿用默认控制 |

点击“运行仿真”后，程序会：

1. 启动或复用 MATLAB Engine，将当前工作目录切换到 `.m` 文件所在位置。
2. 调用 `simulate_arm(struct)` 并传入 GUI 中设定的参数。
3. 在日志中显示仿真状态、采样点数以及最终盘心位置。

### 自定义力表达式示例

- 恒定推力：`0.8`
- 正弦激励：`0.2*sin(2*pi*t)`
- 分段函数（使用 MATLAB 语法）：`(t<2)*0 + (t>=2)*0.5`

> 表达式需要为单变量（时间 `t`），系统会自动构造向量 `[f_1(t); f_2(t); f_3(t)]` 传入 `simulate_arm`。如果某一级需要保持默认策略，可将该级三项全部留空。

若 MATLAB Engine 未安装或启动失败，将在日志与弹窗中给出提示。

## 常见问题

- **MATLAB Engine 导入失败**：请确认已执行 `python setup.py install` 并使用与安装时相同的 Python 解释器。
- **需要新增参数**：可在 `python_ui/app.py` 中扩展 `DEFAULT_PARAMS` 和界面表单，同时更新 `simulate_arm.m` 的参数解析。
- **需要命令行自动化**：可使用 `matlab.engine.start_matlab()` 编写脚本直接调用 `simulate_arm`，无需 GUI。

欢迎根据实际需求继续扩展（例如保存输出数据、加载不同控制策略、调用 Python 绘图库等）。
