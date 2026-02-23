import threading
import tkinter as tk
from tkinter import messagebox, ttk
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

try:
    import matlab.engine  # type: ignore
    import matlab  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    matlab_engine = None
    matlab_module = None
else:  # pragma: no cover - heavy import
    matlab_engine = matlab.engine
    matlab_module = matlab

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PARAMS = {
    "stageCount": 8,
    "duration": 10.0,
    "timeStep": 0.01,
    "derivativeTimeStep": 1e-4,
    "initialHeight": 0.35,
    "enablePlots": False,
}


class MatlabSimulationBridge:
    """轻量封装 MATLAB Engine，负责启动引擎并调用 simulate_arm."""

    def __init__(self, workspace_path: Path) -> None:
        self.workspace_path = workspace_path
        self._engine = None

    def ensure_engine(self) -> None:
        if matlab_engine is None:
            raise RuntimeError(
                (
                    "未检测到 MATLAB Engine for Python。"
                    "请参考 README 安装 matlabengineforpython 并重试。"
                )
            )
        if self._engine is None:
            self._engine = matlab_engine.start_matlab()
            self._engine.cd(str(self.workspace_path), nargout=0)

    def run_simulation(self, params: Dict[str, Any]) -> Any:
        self.ensure_engine()
        matlab_params: Dict[str, Any] = {}
        for key, value in params.items():
            if key == "enablePlots":
                matlab_params[key] = bool(value)
            elif key == "customForceMatrix":
                matlab_params[key] = self._prepare_custom_force(value)
            elif key == "stageCount":
                matlab_params[key] = float(int(value))
            else:
                matlab_params[key] = float(value)
        matlab_params.setdefault("enablePlots", False)
        return self._engine.simulate_arm(matlab_params, nargout=1)

    def _prepare_custom_force(self, value: Any) -> Any:
        if matlab_module is None:
            return value
        if self._is_numeric_nested(value):
            return matlab_module.double(value)
        return value

    @staticmethod
    def _is_numeric_nested(value: Any) -> bool:
        if isinstance(value, (int, float)):
            return True
        if isinstance(value, (list, tuple)) and value:
            return all(
                MatlabSimulationBridge._is_numeric_nested(item)
                for item in value
            )
        if isinstance(value, (list, tuple)):
            return False
        return False


class SimulationApp:
    def __init__(self) -> None:
        self.root = tk.Tk()
        self.root.title("柔性臂仿真控制台")
        self.bridge = MatlabSimulationBridge(WORKSPACE_ROOT)

        self.stage_count_var = tk.IntVar(value=DEFAULT_PARAMS["stageCount"])
        self.duration_var = tk.DoubleVar(value=DEFAULT_PARAMS["duration"])
        self.time_step_var = tk.DoubleVar(value=DEFAULT_PARAMS["timeStep"])
        self.derivative_step_var = tk.DoubleVar(
            value=DEFAULT_PARAMS["derivativeTimeStep"]
        )
        self.initial_height_var = tk.DoubleVar(
            value=DEFAULT_PARAMS["initialHeight"]
        )
        self.enable_plots_var = tk.BooleanVar(
            value=DEFAULT_PARAMS["enablePlots"]
        )
        self.status_var = tk.StringVar(value="准备就绪")

        self.force_entries: List[List[tk.StringVar]] = []
        self.force_refresh_scheduled = False
        self.force_table: Optional[ttk.Frame] = None

        self._build_layout()
        self._refresh_force_inputs()
        self.stage_count_var.trace_add("write", self._on_stage_count_changed)

        if matlab_engine is None:
            self._append_log(
                "未检测到 MATLAB Engine for Python。请先安装后再运行仿真。"
            )

    # --------------------------- UI 构造 --------------------------- #
    def _build_layout(self) -> None:
        main = ttk.Frame(self.root, padding=16)
        main.grid(column=0, row=0, sticky="nsew")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        input_fields = [
            ("级数 (stageCount)", self.stage_count_var),
            ("仿真时长 (s)", self.duration_var),
            ("积分步长 dt (s)", self.time_step_var),
            ("有限差分步长 (s)", self.derivative_step_var),
            ("初始高度 h0 (m)", self.initial_height_var),
        ]

        row = 0
        for label, variable in input_fields:
            row = self._add_labeled_entry(main, row, label, variable)

        plots_btn = ttk.Checkbutton(
            main,
            text="运行时显示 MATLAB 绘图",
            variable=self.enable_plots_var,
        )
        plots_btn.grid(
            column=0,
            row=row,
            columnspan=2,
            sticky="w",
            pady=(8, 0),
        )
        row += 1

        force_group = ttk.LabelFrame(
            main,
            text="各杆驱动力表达式 f_z(t)",
            padding=8,
        )
        force_group.grid(
            column=0,
            row=row,
            columnspan=2,
            sticky="ew",
            pady=(12, 0),
        )
        force_group.columnconfigure(0, weight=1)

        ttk.Label(
            force_group,
            text=(
                "为每级 3 根杆输入沿局部 Z 轴的力函数表达式，"
                "使用变量 t（秒），可调用 MATLAB 内置函数。"
                " 留空表示该级保持默认控制。"
            ),
            wraplength=420,
            justify="left",
        ).grid(column=0, row=0, sticky="w")

        self.force_table = ttk.Frame(force_group)
        self.force_table.grid(column=0, row=1, sticky="ew", pady=(8, 0))
        self.force_table.columnconfigure(0, weight=1)

        row += 1

        ttk.Separator(main).grid(
            column=0,
            row=row,
            columnspan=2,
            sticky="ew",
            pady=12,
        )
        row += 1

        buttons = ttk.Frame(main)
        buttons.grid(column=0, row=row, columnspan=2, sticky="ew")
        buttons.columnconfigure(0, weight=1)
        buttons.columnconfigure(1, weight=1)

        self.run_button = ttk.Button(
            buttons,
            text="运行仿真",
            command=self.run_simulation,
        )
        self.run_button.grid(column=0, row=0, sticky="ew", padx=(0, 4))

        ttk.Button(buttons, text="退出", command=self.root.quit).grid(
            column=1,
            row=0,
            sticky="ew",
            padx=(4, 0),
        )
        row += 1

        ttk.Label(main, textvariable=self.status_var).grid(
            column=0,
            row=row,
            columnspan=2,
            sticky="w",
            pady=(12, 4),
        )
        row += 1

        log_frame = ttk.LabelFrame(main, text="运行日志", padding=8)
        log_frame.grid(column=0, row=row, columnspan=2, sticky="nsew")
        main.rowconfigure(row, weight=1)

        self.log_text = tk.Text(
            log_frame,
            height=12,
            width=60,
            state="disabled",
        )
        self.log_text.pack(fill="both", expand=True)

    def _add_labeled_entry(
        self,
        parent: ttk.Frame,
        row: int,
        label: str,
        variable: tk.Variable,
    ) -> int:
        ttk.Label(parent, text=label).grid(
            column=0,
            row=row,
            sticky="w",
            pady=4,
        )
        entry = ttk.Entry(parent, textvariable=variable)
        entry.grid(column=1, row=row, sticky="ew", pady=4)
        parent.columnconfigure(1, weight=1)
        return row + 1

    def _on_stage_count_changed(self, *_: str) -> None:
        self._schedule_force_refresh()

    def _schedule_force_refresh(self) -> None:
        if self.force_refresh_scheduled:
            return
        self.force_refresh_scheduled = True
        self.root.after(120, self._refresh_force_inputs)

    def _refresh_force_inputs(self) -> None:
        self.force_refresh_scheduled = False
        if self.force_table is None:
            return

        stage_count = max(1, self.stage_count_var.get())
        previous_values = [
            [var.get() for var in stage_vars]
            for stage_vars in self.force_entries
        ]

        for child in self.force_table.winfo_children():
            child.destroy()

        headers = ("级数", "杆1 f_z(t)", "杆2 f_z(t)", "杆3 f_z(t)")
        for col_index, header in enumerate(headers):
            ttk.Label(
                self.force_table,
                text=header,
                anchor="center",
                padding=(4, 2),
                width=16,
            ).grid(
                column=col_index,
                row=0,
                padx=2,
                pady=2,
                sticky="nsew",
            )
            self.force_table.columnconfigure(col_index, weight=1)

        self.force_entries = []
        for stage_index in range(stage_count):
            ttk.Label(
                self.force_table,
                text=f"第{stage_index + 1}级",
                anchor="center",
                padding=(4, 2),
            ).grid(
                column=0,
                row=stage_index + 1,
                padx=2,
                pady=2,
                sticky="nsew",
            )

            row_vars: List[tk.StringVar] = []
            for force_axis in range(3):
                preset = ""
                if stage_index < len(previous_values):
                    stage_values = previous_values[stage_index]
                    if force_axis < len(stage_values):
                        preset = stage_values[force_axis]

                var = tk.StringVar(value=preset)
                entry = ttk.Entry(self.force_table, textvariable=var, width=18)
                entry.grid(
                    column=force_axis + 1,
                    row=stage_index + 1,
                    padx=2,
                    pady=2,
                    sticky="ew",
                )
                row_vars.append(var)

            self.force_entries.append(row_vars)

    def _collect_force_expressions(
        self,
        stage_count: int,
    ) -> Optional[List[str]]:
        if not self.force_entries:
            return None

        expressions: List[str] = []
        has_expression = False

        for stage_index in range(stage_count):
            if stage_index >= len(self.force_entries):
                expressions.append("")
                continue

            row_vars = self.force_entries[stage_index]
            raw_texts = [var.get().strip() for var in row_vars]

            if not any(raw_texts):
                expressions.append("")
                continue

            if any(text == "" for text in raw_texts):
                raise ValueError(
                    f"第{stage_index + 1}级：请同时填写三根杆的力表达式，或全部留空。"
                )

            expression = "[" + "; ".join(raw_texts) + "]"
            expressions.append(expression)
            has_expression = True

        return expressions if has_expression else None

    # --------------------------- 事件逻辑 --------------------------- #
    def run_simulation(self) -> None:
        try:
            params = self._collect_params()
        except ValueError as exc:
            messagebox.showerror("输入错误", str(exc))
            return

        self.status_var.set("正在启动仿真...")
        self.run_button.configure(state=tk.DISABLED)

        thread = threading.Thread(
            target=self._run_in_background,
            args=(params,),
            daemon=True,
        )
        thread.start()

    def _run_in_background(self, params: Dict[str, Any]) -> None:
        try:
            result = self.bridge.run_simulation(params)
        except Exception as exc:  # noqa: BLE001 - 向 UI 回传
            self.root.after(0, self._handle_failure, exc)
        else:
            self.root.after(0, self._handle_success, params, result)

    def _handle_success(self, params: Dict[str, Any], result: Any) -> None:
        stage_count = params["stageCount"]
        duration = params["duration"]
        final_position = self._extract_final_position(result)
        trace_len = self._extract_trace_length(result)

        summary = (
            f"仿真完成：级数 {stage_count}，时长 {duration:g} s，"
            f"采样点 {trace_len}，最终盘心位置 {final_position}"
        )
        self._append_log(summary)
        self.status_var.set("仿真完成")
        self.run_button.configure(state=tk.NORMAL)

    def _handle_failure(self, exc: Exception) -> None:
        self.status_var.set("运行失败")
        self.run_button.configure(state=tk.NORMAL)
        messagebox.showerror("运行失败", str(exc))
        self._append_log(f"错误：{exc}")

    # --------------------------- 参数与结果处理 --------------------------- #
    def _collect_params(self) -> Dict[str, Any]:
        stage_count = self.stage_count_var.get()
        if stage_count <= 0:
            raise ValueError("级数必须为正整数")

        self._refresh_force_inputs()

        duration = self._ensure_positive(self.duration_var.get(), "仿真时长")
        time_step = self._ensure_positive(self.time_step_var.get(), "积分步长")
        derivative_step = self._ensure_positive(
            self.derivative_step_var.get(),
            "有限差分步长",
        )

        params: Dict[str, Any] = {
            "stageCount": stage_count,
            "duration": duration,
            "timeStep": time_step,
            "derivativeTimeStep": derivative_step,
            "initialHeight": self.initial_height_var.get(),
            "enablePlots": self.enable_plots_var.get(),
        }

        custom_force = self._collect_force_expressions(stage_count)
        if custom_force is not None:
            params["customForceMatrix"] = custom_force

        return params

    def _ensure_positive(self, value: float, label: str) -> float:
        if value <= 0:
            raise ValueError(f"{label} 必须为正值")
        return value

    def _extract_final_position(self, result: Any) -> List[float]:
        matrix = self._get_struct_field(result, "finalPosition")
        if not matrix:
            return []
        last_row = matrix[-1]
        return [round(float(component), 6) for component in last_row]

    def _extract_trace_length(self, result: Any) -> int:
        time_series = self._get_struct_field(result, "time")
        return len(time_series)

    def _get_struct_field(self, result: Any, field_name: str) -> List[Any]:
        try:
            field = result[field_name]
        except (TypeError, KeyError):
            return []
        return self._to_python_numeric(field)

    def _to_python_numeric(self, value: Any) -> Any:
        if isinstance(value, (float, int, bool)):
            return value
        if hasattr(value, "tolist"):
            converted = value.tolist()
        else:
            converted = value
        if isinstance(converted, Iterable) and not isinstance(
            converted,
            (str, bytes),
        ):
            return [self._to_python_numeric(item) for item in converted]
        return converted

    def _append_log(self, message: str) -> None:
        self.log_text.configure(state="normal")
        self.log_text.insert("end", message + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    # --------------------------- 入口 --------------------------- #
    def run(self) -> None:
        self.root.mainloop()


def main() -> None:
    app = SimulationApp()
    app.run()


if __name__ == "__main__":
    main()
