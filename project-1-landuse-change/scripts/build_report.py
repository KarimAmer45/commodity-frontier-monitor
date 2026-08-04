"""Rebuild the portfolio report without loading geospatial analysis packages."""

import ast
import json
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PDF_DIR = ROOT / "output" / "pdf"


def load_report_builder():
    """Load only the report functions from the full analysis script."""
    source_path = Path(__file__).with_name("run_local_analysis.py")
    tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
    functions = [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name in {"footer", "build_report_pdf"}
    ]
    module = ast.Module(body=functions, type_ignores=[])
    exec(compile(module, str(source_path), "exec"), globals())


def main():
    OUTPUT_PDF_DIR.mkdir(parents=True, exist_ok=True)
    load_report_builder()
    payload = json.loads((ROOT / "results.json").read_text(encoding="utf-8"))
    results = payload.get("results", payload)
    output = build_report_pdf(results, ROOT / "map" / "landuse_change_map.png")
    print(output)


if __name__ == "__main__":
    main()
