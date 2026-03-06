#!/usr/bin/env python3
import argparse
import contextlib
import html
import io
import json
import os
import traceback
from pathlib import Path


def _as_text(source):
    if isinstance(source, list):
        return "".join(source)
    return str(source)


def _execute_with_nbconvert(
    notebook_path: Path,
    run_dir: Path,
    out_notebook_path: Path,
    out_html_path: Path,
    timeout_s: int,
):
    import nbformat  # type: ignore
    from nbconvert import HTMLExporter  # type: ignore
    from nbconvert.preprocessors import ExecutePreprocessor  # type: ignore

    nb = nbformat.read(str(notebook_path), as_version=4)

    # Notebook derives RUN_DIR from Path.cwd(); execute from run_dir.
    ep = ExecutePreprocessor(timeout=timeout_s, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": str(run_dir)}})

    with out_notebook_path.open("w", encoding="utf-8") as f:
        nbformat.write(nb, f)

    exporter = HTMLExporter()
    body, _ = exporter.from_notebook_node(nb)
    out_html_path.write_text(body, encoding="utf-8")


def _fallback_execute_notebook(
    notebook_path: Path,
    run_dir: Path,
    out_notebook_path: Path,
    out_html_path: Path,
):
    nb = json.loads(notebook_path.read_text(encoding="utf-8"))
    cells = nb.get("cells", [])

    old_cwd = Path.cwd()
    os.chdir(run_dir)

    env = {"__name__": "__main__"}

    def display(obj):
        print(obj)

    env["display"] = display

    try:
        exec(
            "import matplotlib\n"
            "matplotlib.use('Agg')\n",
            env,
            env,
        )
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Fallback execution requires matplotlib/pandas/numpy in the active Python environment."
        ) from exc

    exec_count = 1
    try:
        for idx, cell in enumerate(cells):
            if cell.get("cell_type") != "code":
                continue

            source = _as_text(cell.get("source", ""))
            stdout = io.StringIO()
            stderr = io.StringIO()
            outputs = []

            try:
                with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                    exec(source, env, env)
            except Exception as exc:
                std_out_text = stdout.getvalue()
                std_err_text = stderr.getvalue()
                if std_out_text:
                    outputs.append({"name": "stdout", "output_type": "stream", "text": std_out_text})
                if std_err_text:
                    outputs.append({"name": "stderr", "output_type": "stream", "text": std_err_text})

                tb_lines = traceback.format_exc().splitlines()
                outputs.append(
                    {
                        "ename": type(exc).__name__,
                        "evalue": str(exc),
                        "output_type": "error",
                        "traceback": tb_lines,
                    }
                )
                cell["execution_count"] = exec_count
                cell["outputs"] = outputs
                raise RuntimeError(f"Notebook execution failed at code cell index {idx}") from exc

            std_out_text = stdout.getvalue()
            std_err_text = stderr.getvalue()
            if std_out_text:
                outputs.append({"name": "stdout", "output_type": "stream", "text": std_out_text})
            if std_err_text:
                outputs.append({"name": "stderr", "output_type": "stream", "text": std_err_text})

            cell["execution_count"] = exec_count
            cell["outputs"] = outputs
            exec_count += 1
    finally:
        os.chdir(old_cwd)

    out_notebook_path.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    out_html_path.write_text(_fallback_html(nb, run_dir), encoding="utf-8")


def _fallback_html(nb, run_dir: Path):
    report_dir = run_dir / "report"
    fig_dir = report_dir / "figures"

    parts = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'>",
        "<title>simio single run report</title>",
        "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;}pre{background:#f4f4f4;padding:0.75rem;overflow-x:auto;}h1,h2,h3{margin-top:1.5rem;} .cell{margin:1rem 0;} .out{margin-top:0.5rem;}</style>",
        "</head><body>",
        f"<h1>simio single run report</h1><p>Run directory: <code>{html.escape(str(run_dir))}</code></p>",
    ]

    for cell in nb.get("cells", []):
        ctype = cell.get("cell_type")
        if ctype == "markdown":
            md = _as_text(cell.get("source", ""))
            parts.append(f"<div class='cell'><pre>{html.escape(md)}</pre></div>")
            continue
        if ctype != "code":
            continue

        src = _as_text(cell.get("source", ""))
        parts.append("<div class='cell'>")
        parts.append(f"<pre><code>{html.escape(src)}</code></pre>")
        for out in cell.get("outputs", []):
            otype = out.get("output_type")
            if otype == "stream":
                txt = out.get("text", "")
                if isinstance(txt, list):
                    txt = "".join(txt)
                parts.append(f"<pre class='out'>{html.escape(str(txt))}</pre>")
            elif otype == "error":
                tb = out.get("traceback", [])
                parts.append(f"<pre class='out'>{html.escape(chr(10).join(tb))}</pre>")
        parts.append("</div>")

    if fig_dir.exists():
        pngs = sorted(fig_dir.glob("*.png"))
        if pngs:
            parts.append("<h2>Figures</h2>")
            for png in pngs:
                rel = png.relative_to(report_dir)
                parts.append(f"<h3>{html.escape(png.name)}</h3>")
                parts.append(f"<img src='{html.escape(str(rel))}' style='max-width:100%;border:1px solid #ddd;padding:4px;'>")

    summary = report_dir / "summary_table.csv"
    if summary.exists():
        parts.append("<h2>Summary Table</h2>")
        parts.append(f"<p><a href='{html.escape(summary.name)}'>{html.escape(summary.name)}</a></p>")

    parts.append("</body></html>")
    return "\n".join(parts)


def main():
    parser = argparse.ArgumentParser(
        description="Execute notebooks/simio_single_run_report.ipynb for one run directory and export notebook + HTML."
    )
    parser.add_argument("run_dir", help="Run directory containing simio CSV outputs")
    parser.add_argument(
        "--notebook",
        default="notebooks/simio_single_run_report.ipynb",
        help="Template notebook path",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=1800,
        help="Execution timeout in seconds (nbconvert path only)",
    )
    parser.add_argument(
        "--force-fallback",
        action="store_true",
        help="Use stdlib fallback executor even if nbconvert stack is installed",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    notebook_path = (repo_root / args.notebook).resolve()
    run_dir = Path(args.run_dir).resolve()
    report_dir = run_dir / "report"
    fig_dir = report_dir / "figures"

    if not notebook_path.exists():
        raise FileNotFoundError(f"Notebook template not found: {notebook_path}")
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    report_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    out_notebook_path = report_dir / "simio_single_run_report.ipynb"
    out_html_path = report_dir / "simio_single_run_report.html"

    if args.force_fallback:
        _fallback_execute_notebook(notebook_path, run_dir, out_notebook_path, out_html_path)
        print("Rendered report with fallback executor.")
    else:
        try:
            _execute_with_nbconvert(
                notebook_path=notebook_path,
                run_dir=run_dir,
                out_notebook_path=out_notebook_path,
                out_html_path=out_html_path,
                timeout_s=args.timeout,
            )
            print("Rendered report with nbconvert/nbclient stack.")
        except (ModuleNotFoundError, PermissionError):
            print("nbconvert execution unavailable in this environment; using fallback executor.")
            _fallback_execute_notebook(notebook_path, run_dir, out_notebook_path, out_html_path)

    print(f"Notebook: {out_notebook_path}")
    print(f"HTML    : {out_html_path}")
    print(f"Figures : {fig_dir}")
    print(f"Summary : {report_dir / 'summary_table.csv'}")


if __name__ == "__main__":
    main()
