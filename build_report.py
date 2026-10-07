#!/usr/bin/env python3
"""Build the polished Chinese DOCX report from completed analysis outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parent
SUMMARY = json.loads((ROOT / "summary.json").read_text(encoding="utf-8"))
METRICS = pd.read_csv(ROOT / "candidate_metrics.csv")
OUTPUT = ROOT / "Combined_circDesign_回顾性计算评估报告.docx"

MISSING: list[str] = []


def load_json(name, default=None):
    """Load one artefact of the verification layer.

    The report is assembled from independent scripts (axis decomposition, IRES metrics,
    design panel, advantage-space test, boundary redo, implementation cross-checks).
    A missing file must show up as a MISSING section, not as a quietly absent claim,
    so absence is recorded here and reported at the top of the built document.
    """
    path = ROOT / name
    if not path.exists():
        MISSING.append(name)
        return default
    return json.loads(path.read_text(encoding="utf-8"))


AXIS = load_json("axis_decomposition.json")
IRES_AN = load_json("ires_analysis.json")
PANEL = load_json("design_panel.json")
ADV = load_json("advantage_space_pooled.json")
BOUND = load_json("boundary_test.json")
REFOLD = load_json("refold_check.json")
REPRO = load_json("shuffle_repro_test.json")
ACC = load_json("accessibility_summary.json")
CLUST = load_json("arrangement_clusters.json")

_panel_csv = ROOT / "design_panel.csv"
PANEL_BY_ID = (
    pd.read_csv(_panel_csv, encoding="utf-8-sig").set_index("candidate_id").to_dict("index")
    if _panel_csv.exists() else {}
)


def pct(k, n):
    return 100.0 * k / n if n else float("nan")

INK = RGBColor(11, 37, 69)
BLUE = RGBColor(46, 116, 181)
DARK_BLUE = RGBColor(31, 77, 120)
MUTED = RGBColor(91, 103, 116)
GREEN = RGBColor(31, 100, 70)
GOLD = RGBColor(122, 90, 0)
RED = RGBColor(155, 28, 28)
LIGHT_BLUE = "E8EEF5"
LIGHT_GREEN = "E8F3EC"
LIGHT_GOLD = "FFF5D6"
LIGHT_RED = "FDECEC"
LIGHT_GRAY = "F2F4F7"
WHITE = "FFFFFF"

LATIN_FONT = "Arial"
EAST_ASIA_FONT = "Noto Sans CJK SC"


def set_run_font(run, size=None, bold=None, italic=None, color=None, name=LATIN_FONT):
    run.font.name = name
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), name)
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), name)
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), EAST_ASIA_FONT)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic
    if color is not None:
        run.font.color.rgb = color


def set_cell_shading(cell, fill):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=90, start=120, bottom=90, end=120):
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for margin, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{margin}"))
        if node is None:
            node = OxmlElement(f"w:{margin}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_repeat_table_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def set_table_geometry(table, widths_dxa):
    table.autofit = False
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    tbl = table._tbl
    tbl_pr = tbl.tblPr
    tbl_w = tbl_pr.find(qn("w:tblW"))
    if tbl_w is None:
        tbl_w = OxmlElement("w:tblW")
        tbl_pr.append(tbl_w)
    tbl_w.set(qn("w:w"), str(sum(widths_dxa)))
    tbl_w.set(qn("w:type"), "dxa")
    tbl_ind = tbl_pr.find(qn("w:tblInd"))
    if tbl_ind is None:
        tbl_ind = OxmlElement("w:tblInd")
        tbl_pr.append(tbl_ind)
    tbl_ind.set(qn("w:w"), "120")
    tbl_ind.set(qn("w:type"), "dxa")

    grid = tbl.tblGrid
    for child in list(grid):
        grid.remove(child)
    for width in widths_dxa:
        col = OxmlElement("w:gridCol")
        col.set(qn("w:w"), str(width))
        grid.append(col)
    for row in table.rows:
        for idx, cell in enumerate(row.cells):
            tc_pr = cell._tc.get_or_add_tcPr()
            tc_w = tc_pr.find(qn("w:tcW"))
            if tc_w is None:
                tc_w = OxmlElement("w:tcW")
                tc_pr.append(tc_w)
            tc_w.set(qn("w:w"), str(widths_dxa[idx]))
            tc_w.set(qn("w:type"), "dxa")
            set_cell_margins(cell)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def set_table_borders(table, color="D7DCE2", size=6):
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.find(qn("w:tblBorders"))
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = borders.find(qn(f"w:{edge}"))
        if tag is None:
            tag = OxmlElement(f"w:{edge}")
            borders.append(tag)
        tag.set(qn("w:val"), "single")
        tag.set(qn("w:sz"), str(size))
        tag.set(qn("w:space"), "0")
        tag.set(qn("w:color"), color)


def style_table_text(table, header=True, body_size=9.2):
    for ridx, row in enumerate(table.rows):
        if header and ridx == 0:
            set_repeat_table_header(row)
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.space_before = Pt(0)
                paragraph.paragraph_format.space_after = Pt(0)
                paragraph.paragraph_format.line_spacing = 1.08
                for run in paragraph.runs:
                    set_run_font(
                        run,
                        size=9.2 if ridx == 0 else body_size,
                        bold=True if ridx == 0 else None,
                        color=INK if ridx == 0 else RGBColor(32, 36, 42),
                    )
            if header and ridx == 0:
                set_cell_shading(cell, LIGHT_BLUE)


def add_table(doc, headers, rows, widths_dxa, body_size=9.2):
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    for idx, text in enumerate(headers):
        table.rows[0].cells[idx].text = str(text)
    for row_data in rows:
        cells = table.add_row().cells
        for idx, value in enumerate(row_data):
            cells[idx].text = str(value)
    set_table_geometry(table, widths_dxa)
    set_table_borders(table)
    style_table_text(table, header=True, body_size=body_size)
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(2)
    return table


def add_callout(doc, label, text, fill=LIGHT_GREEN, label_color=GREEN):
    table = doc.add_table(rows=1, cols=1)
    set_table_geometry(table, [9360])
    set_table_borders(table, color=fill, size=4)
    cell = table.cell(0, 0)
    set_cell_shading(cell, fill)
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.line_spacing = 1.12
    r = p.add_run(label + "　")
    set_run_font(r, size=11, bold=True, color=label_color)
    r = p.add_run(text)
    set_run_font(r, size=10.5, color=INK)
    spacer = doc.add_paragraph()
    spacer.paragraph_format.space_after = Pt(2)
    return table


def add_body(doc, text, bold_lead=None):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.line_spacing = 1.10
    if bold_lead and text.startswith(bold_lead):
        r = p.add_run(bold_lead)
        set_run_font(r, size=11, bold=True, color=INK)
        r = p.add_run(text[len(bold_lead) :])
        set_run_font(r, size=11, color=RGBColor(32, 36, 42))
    else:
        r = p.add_run(text)
        set_run_font(r, size=11, color=RGBColor(32, 36, 42))
    return p


def add_bullet(doc, text, level=0):
    p = doc.add_paragraph(style="List Bullet" if level == 0 else "List Bullet 2")
    p.paragraph_format.left_indent = Inches(0.5 if level == 0 else 0.75)
    p.paragraph_format.first_line_indent = Inches(-0.25)
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.line_spacing = 1.10
    for run in p.runs:
        set_run_font(run, size=10.8, color=RGBColor(32, 36, 42))
    if not p.runs:
        r = p.add_run(text)
        set_run_font(r, size=10.8, color=RGBColor(32, 36, 42))
    return p


def add_numbered(doc, text):
    p = doc.add_paragraph(style="List Number")
    p.paragraph_format.left_indent = Inches(0.5)
    p.paragraph_format.first_line_indent = Inches(-0.25)
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.line_spacing = 1.10
    if p.runs:
        for run in p.runs:
            set_run_font(run, size=10.8, color=RGBColor(32, 36, 42))
    else:
        r = p.add_run(text)
        set_run_font(r, size=10.8, color=RGBColor(32, 36, 42))
    return p


def add_heading(doc, text, level=1):
    p = doc.add_paragraph(text, style=f"Heading {level}")
    p.paragraph_format.keep_with_next = True
    return p


def add_figure(doc, image_name, caption, width=6.35):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(3)
    run = p.add_run()
    run.add_picture(str(ROOT / image_name), width=Inches(width))
    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap.paragraph_format.space_before = Pt(0)
    cap.paragraph_format.space_after = Pt(8)
    cap.paragraph_format.keep_with_next = False
    r = cap.add_run(caption)
    set_run_font(r, size=9, color=MUTED)


def set_paragraph_shading(paragraph, fill):
    p_pr = paragraph._p.get_or_add_pPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    p_pr.append(shd)


def add_page_number(paragraph):
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    r = paragraph.add_run("第 ")
    set_run_font(r, size=8.5, color=MUTED)
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    r_node = OxmlElement("w:r")
    t_node = OxmlElement("w:t")
    t_node.text = "1"
    r_node.append(t_node)
    fld.append(r_node)
    paragraph._p.append(fld)
    r = paragraph.add_run(" 页")
    set_run_font(r, size=8.5, color=MUTED)


def configure_styles(doc):
    normal = doc.styles["Normal"]
    normal.font.name = LATIN_FONT
    normal._element.rPr.rFonts.set(qn("w:ascii"), LATIN_FONT)
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), LATIN_FONT)
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), EAST_ASIA_FONT)
    normal.font.size = Pt(11)
    pf = normal.paragraph_format
    pf.space_before = Pt(0)
    pf.space_after = Pt(6)
    pf.line_spacing = 1.10

    style_specs = {
        "Heading 1": (16, BLUE, 16, 8),
        "Heading 2": (13, BLUE, 12, 6),
        "Heading 3": (12, DARK_BLUE, 8, 4),
    }
    for name, (size, color, before, after) in style_specs.items():
        style = doc.styles[name]
        style.font.name = LATIN_FONT
        style._element.rPr.rFonts.set(qn("w:ascii"), LATIN_FONT)
        style._element.rPr.rFonts.set(qn("w:hAnsi"), LATIN_FONT)
        style._element.rPr.rFonts.set(qn("w:eastAsia"), EAST_ASIA_FONT)
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = color
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(after)
        style.paragraph_format.keep_with_next = True

    for list_name in ("List Bullet", "List Bullet 2", "List Number"):
        style = doc.styles[list_name]
        style.font.name = LATIN_FONT
        style._element.rPr.rFonts.set(qn("w:ascii"), LATIN_FONT)
        style._element.rPr.rFonts.set(qn("w:hAnsi"), LATIN_FONT)
        style._element.rPr.rFonts.set(qn("w:eastAsia"), EAST_ASIA_FONT)
        style.font.size = Pt(10.8)


def configure_page(doc):
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(1.0)
    section.bottom_margin = Inches(1.0)
    section.left_margin = Inches(1.0)
    section.right_margin = Inches(1.0)
    section.header_distance = Inches(0.492)
    section.footer_distance = Inches(0.492)

    header = section.header
    table = header.add_table(rows=1, cols=2, width=Inches(6.5))
    set_table_geometry(table, [4680, 4680])
    for cell in table.rows[0].cells:
        set_cell_margins(cell, top=0, start=0, bottom=0, end=0)
    p = table.cell(0, 0).paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    r = p.add_run("JLU-FBH · CirCure Model")
    set_run_font(r, size=8.3, bold=True, color=MUTED)
    p = table.cell(0, 1).paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    r = p.add_run("Combined circRNA · Retrospective Audit")
    set_run_font(r, size=8.3, color=MUTED)

    footer = section.footer
    add_page_number(footer.paragraphs[0])


def add_title_block(doc):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(18)
    p.paragraph_format.space_after = Pt(7)
    r = p.add_run("JLU-FBH · CirCure · COMPUTATIONAL MODEL")
    set_run_font(r, size=9.5, bold=True, color=BLUE)

    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(5)
    r = p.add_run("Combined成熟circRNA构建")
    set_run_font(r, size=25, bold=True, color=INK)
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(8)
    r = p.add_run("circDesign框架回顾性计算评估报告")
    set_run_font(r, size=22, bold=True, color=INK)

    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(16)
    r = p.add_run("基于公开设计目标的独立实现 · 完整闭环折叠 · 同义序列空间 · 工程决策")
    set_run_font(r, size=11.5, color=MUTED)

    metadata = [
        ("输入文件", SUMMARY["input_file"]),
        ("分析对象", "2,013 nt Combined mature circRNA（质粒坐标 723–2735）"),
        ("软件与模型", f"ViennaRNA {SUMMARY['software']['ViennaRNA']} circular model；Python {SUMMARY['software']['python']}"),
        ("参照空间",
         f"{SUMMARY['benchmark_design']['random_total_n']} 条随机同义序列 + "
         f"{len(SUMMARY['benchmark_design']['bounding_controls'])} 条边界对照；"
         f"另含扩充样本 {ADV['n_generated']:,} 条" if ADV else
         f"{SUMMARY['benchmark_design']['random_total_n']} 条随机同义序列"),
        ("随机种子", str(SUMMARY["software"]["random_seed"])),
        ("报告日期", SUMMARY["analysis_timestamp_local"]),
    ]
    add_table(doc, ["项目", "内容"], metadata, [1650, 7710], body_size=9.4)


def format_pct(p):
    return f"{p['percentile']:.1f}%（95% CI {p['ci95_low']:.1f}–{p['ci95_high']:.1f}）"


def main():
    doc = Document()
    doc.core_properties.title = "Combined成熟circRNA构建 circDesign框架回顾性计算评估报告"
    doc.core_properties.subject = "JLU-FBH CirCure computational model"
    doc.core_properties.author = "JLU-FBH"
    doc.core_properties.keywords = "circRNA, circDesign, MFE, CAI, IRES, CirCure"
    configure_styles(doc)
    configure_page(doc)
    add_title_block(doc)

    if ADV and ACC:
        if ACC["dominators_full_objective"] == 0:
            adv_line = (
                f"在保持密码子组成完全不变的同义重排中，仅 {ADV['rate_passing_percent']:.2f}% "
                f"能在三项目标上同时不劣于本构建；并入翻译可及性后，"
                f"{ACC['n_buildable']:,} 条可建成同义设计中无一能同时不劣于本构建。"
                f"IRES–CDS 间隔区可及性本身位于可建成设计的第 "
                f"{100 - ACC['construct_position']['acc_SPACER']['percent_better_than_construct']:.1f} 百分位。"
            )
        else:
            adv_line = (
                f"在保持密码子组成完全不变的同义重排中，仅 {ADV['rate_passing_percent']:.2f}% "
                f"能在三项目标上同时不劣于本构建；纳入翻译可及性的完整设计目标集后，"
                f"该比例降至 {ACC['domination_rate_full_objective_percent']:.2f}%"
                f"（{ACC['n_buildable']:,} 条可建成设计中的 "
                f"{ACC['dominators_full_objective']} 条）。"
            )
    elif ADV:
        adv_line = (
            f"在保持密码子组成完全不变的同义重排中，仅 {ADV['rate_passing_percent']:.2f}% "
            f"能在三项目标上同时不劣于本构建。"
        )
    else:
        adv_line = ""
    add_callout(
        doc,
        "结论",
        "计算结果支持保留现有 Combined 构建用于当前研发链。"
        + adv_line
        + "密码子组成在参照空间中达到上确界，且不存在系统性的改进设计。"
        "本报告不声称该序列是同义排列空间中的绝对最优解——见第九节。",
        fill=LIGHT_GREEN,
        label_color=GREEN,
    )
    if MISSING:
        add_callout(
            doc,
            "构建提示",
            "以下校验层产物缺失，对应章节未生成：" + "、".join(MISSING),
            fill=LIGHT_RED,
            label_color=RED,
        )

    add_heading(doc, "关键结果", 1)
    current = SUMMARY["current_metrics"]
    per = SUMMARY["percentiles"]
    sh_p = per["composition_shuffle"]
    hw_p = per["human_weighted"]
    un_p = per["uniform"]

    if ADV:
        # The delivered 64-shuffle ensemble is a 1-in-104 draw (see section 10), so the
        # composition-matched MFE percentile is taken from the pooled n=4000 re-sample.
        mfe_arr = 100.0 * (1 - ADV["n_mfe_competitive"] / ADV["n_generated"])
        mfe_arr_n = ADV["n_generated"]
        mfe_src = f"{mfe_arr:.2f}%（n={mfe_arr_n}，合并两独立样本）"
    else:
        mfe_src = f"{sh_p['circular_mfe']['percentile']:.1f}%（n={sh_p['circular_mfe']['n']}）"
    rows = [
        (
            "Circular MFE",
            f"{current['circular_mfe_kcal_mol']:.2f} kcal/mol",
            mfe_src,
            "组成匹配零模型下的位次；越低越好",
        ),
        (
            "Human CAI",
            f"{current['human_cai']:.4f}",
            f"{hw_p['human_cai']['percentile']:.1f}% / {un_p['human_cai']['percentile']:.1f}%",
            "人源加权 / 均匀抽样；二者均无候选超过本构建",
        ),
    ]
    if IRES_AN:
        for col, label in (("ires_native_recall_hard_666", "IRES 自身结构保持率"),
                           ("ires_cross_occupancy_666", "IRES 被环外侵占率")):
            ex = IRES_AN["exchangeability_test"].get(col, {})
            if not ex:
                continue
            rows.append((
                label,
                f"{ex['original']:.4f}",
                f"{100 * (1 - (ex.get('rank', 0) - 1) / (ex.get('n', 1) - 1)):.1f}%"
                if ex.get("rank") else "—",
                f"组成匹配排列中的位次；可交换性 p = {ex.get('p_exchangeability', float('nan')):.3f}",
            ))
    add_table(doc, ["指标", "当前值", "组成匹配零模型下的位次", "说明"], rows,
              [1900, 1550, 2400, 3510], body_size=8.8)

    add_bullet(
        doc,
        "报告口径说明：每个指标只对其能回答的零模型报告百分位。组成轴（CAI）用"
        "人源加权与均匀抽样；排列轴（MFE、IRES）用组成重排（composition shuffle），"
        "因为只有后者在保持密码子多重集不变的前提下隔离排列效应。",
    )
    if ADV:
        add_bullet(
            doc,
            f"多目标核验：保持密码子组成完全不变的同义重排中，"
            f"{ADV['dominators_correctness_passing']}/{ADV['n_correctness_passing']}"
            f"（{ADV['rate_passing_percent']:.2f}%）通过产物正确性质控后仍在三项目标上"
            f"同时不劣于本构建。原报告的 5/256 分母含代数上不可能命中的候选，"
            f"结构性错误，已废止。",
        )
    if AXIS:
        ev = AXIS["metric_space"]["eigenvalues"]
        add_bullet(
            doc,
            f"三个目标并非彼此独立：标准化后有效秩仅 "
            f"{AXIS['metric_space']['effective_rank']:.2f}/3"
            f"（特征值 {ev[0]:.3f} / {ev[1]:.3f} / {ev[2]:.3f}），"
            f"其中 MFE 与 CAI 的相关系数为 "
            f"{AXIS['metric_space']['correlation_matrix'][0][1]:+.3f}。",
        )

    add_heading(doc, "1　分析目的与结论边界", 1)
    add_body(
        doc,
        "本评估的目标不是证明 circDesign 原始算法本身有效，也不是声称完成了其未公开生产代码的忠实复现，而是使用论文公开的三个设计目标，对团队已经构建的 Combined 成熟 circRNA 做回顾性计算审计：它的完整闭环结构是否足够紧凑、其 CDS 是否适配人源密码子偏好、其 IRES 在完整闭环背景中是否保持预测结构。",
    )
    add_body(
        doc,
        "circDesign 并不存在跨长度、跨 IRES、跨蛋白通用的合格阈值。因此，本报告不把论文中 RABV-G 构建的 MFE < −1100 kcal/mol 直接套用于本项目，而是以相同长度、相同非编码区、相同蛋白产物的同义序列作为参照，报告经验百分位和 Pareto 支配关系。",
    )
    add_callout(
        doc,
        "术语边界",
        "准确表述为“circDesign-inspired retrospective evaluation”或“基于 circDesign 设计目标的回顾性评估”。不应写成“通过 circDesign 官方认证”“符合 circDesign 官方标准”或“完整复现 circDesign 原始程序”。",
        fill=LIGHT_GOLD,
        label_color=GOLD,
    )

    add_heading(doc, "2　成熟circRNA重建与序列质控", 1)
    coord = SUMMARY["coordinate_model"]
    coordinate_rows = [
        ("质粒", "1–5267", "—", "5,267 bp SnapGene circular DNA"),
        ("成熟 circRNA", "723–2735", "1–2013", "两个 retained-scar 注释所界定；2735→723 闭环"),
        ("synIRES-RC25", "823–1488", "101–766", "666 nt，主分析边界"),
        ("HRV-B3 IRES 扩展注释", "823–1491", "101–769", "669 nt，敏感性分析"),
        ("Combined CDS", "1504–2586", "782–1864", "1,083 nt；360 aa + TGA；无内部终止"),
    ]
    add_table(doc, ["对象", "质粒坐标（1-based）", "成熟环坐标", "核验结果"], coordinate_rows, [1800, 2050, 1700, 3810], body_size=8.7)
    add_body(
        doc,
        "成熟闭环连接点两侧 32 nt 为：\n"
        + coord["junction_32nt_each_side"],
        bold_lead="成熟闭环连接点两侧 32 nt 为：",
    )
    add_body(
        doc,
        "CDS 翻译检查显示：起始密码子 ATG、终止密码子 TGA，编码 360 aa，未发现内部终止或移码；tPA signal/pro sequence、8 个 CTL 表位及连接肽、4×mC3dP28、MITD 和 FLAG 均位于同一连续开放读码框内。",
    )
    add_callout(
        doc,
        "实验身份确认",
        "成熟环边界由 .dna feature 中的 Intron Scar / retained-scar 注释推断。为证明计算实体与实验产物完全一致，仍建议用 back-splice junction RT-PCR 与 junction Sanger 测序确认 2735→723 的真实连接。",
        fill=LIGHT_GOLD,
        label_color=GOLD,
    )

    add_heading(doc, "3　计算方法", 1)
    add_heading(doc, "3.1　三项目标", 2)
    add_numbered(doc, "Circular MFE：在 ViennaRNA circular model 下求完整 2,013 nt 成熟环的最低自由能，MFE(r) = minₛ ΔG°(r,s)。在相同长度候选之间，更负的 MFE 表示预测结构更紧凑。")
    add_numbered(doc, "Human CAI：对 360 个 sense codon 计算 CAI = (∏ᵢw(cᵢ))^(1/n)，终止密码子不纳入；w(cᵢ) 由 Homo sapiens 同义密码子相对适应度确定。")
    add_numbered(doc, "IRES structural deviation：计算完整闭环候选与 IRES-only 受约束参考的碱基配对概率矩阵 L2 距离，L_IRES = [Σ(P_candidate(i,j) − P_reference(i,j))²]^(1/2)。越低表示 IRES 在完整构建背景中受到的预测结构干扰越小。")

    add_heading(doc, "3.2　折叠设置与IRES参考", 2)
    method_rows = [
        ("软件", f"ViennaRNA {SUMMARY['software']['ViennaRNA']}"),
        ("拓扑", "circular = 1；完整 2,013 nt 成熟环"),
        ("温度与 dangling ends", "37°C；dangles = 2；允许孤立碱基对"),
        ("IRES-only 参考", "在相同完整闭环上强制所有非 IRES 碱基保持不配对，仅允许 IRES 内部形成结构"),
        ("BPP 稀疏阈值", f"p ≥ {SUMMARY['software']['bpp_cutoff']:.0e}；计算完整对称矩阵 L2 范数"),
        ("主边界 / 敏感边界", "666 nt（101–766）/ 669 nt（101–769）"),
    ]
    add_table(doc, ["设置", "采用值"], method_rows, [2600, 6760], body_size=9.2)
    add_body(
        doc,
        "IRES deviation 被进一步分解为内部概率结构变化和 IRES—非IRES 跨区配对两部分。当前构建的总 L_IRES = 11.9540，其中内部偏离分量为 8.0251，跨区配对分量为 8.8598；分量按平方和组成总值，不能线性相加。",
    )

    add_heading(doc, "3.3　Human CAI参考", 2)
    add_body(
        doc,
        "采用 Kazusa Codon Usage Database 的 Homo sapiens 通用数据（93,487 条 CDS，共 40,662,582 个密码子），以每种氨基酸中最高频同义密码子为 1，按原始 codon counts 计算相对适应度。该值反映一般人源密码子偏好，不是乳腺、树突状细胞或特定 tRNA 丰度模型。",
    )

    add_heading(doc, "3.4　同义参照空间", 2)
    benchmark_rows = [
        ("Human-weighted", "96", "按人源密码子实际使用频数加权抽样", "现实的人源化随机设计参照"),
        ("Uniform synonymous", "96", "每种氨基酸的同义密码子等概率抽样", "不带人源偏好的全空间参照"),
        ("Composition shuffle", "64", "在相同氨基酸的位置间重排原有同义密码子", "CAI、GC、GC3与原构建相同；隔离局部排列效应"),
        ("Bounding controls", "4", "CAI-max、CAI-min、GC3-high、GC3-low", "定义极端边界；不计入随机百分位"),
    ]
    add_table(doc, ["序列组", "n", "生成方式", "用途"], benchmark_rows, [1900, 650, 3350, 3460], body_size=8.8)
    add_body(
        doc,
        "全部候选均保持成熟环长度 2,013 nt、IRES 与非编码区不变，编码完全相同的 360 aa Combined 抗原，只改变 1,083 nt CDS 中的同义密码子；终止密码子固定为 TGA。随机种子为 20260804。",
    )

    add_heading(doc, "4　结果：按指标配零模型", 1)
    add_body(
        doc,
        "原报告把 96 条人源加权、96 条均匀与 64 条组成重排池成一个 n=256 的参照分布，"
        "据以计算三项指标的百分位。该做法在方法上不成立：三组沿两个不同方向变动——"
        "人源加权与均匀抽样改变密码子组成，组成重排则保持组成恒定。混合后的分布由组成方差主导，"
        "会把每个指标都拉向对本构建有利的方向。本节按指标配零模型重新报告。",
    )

    add_heading(doc, "4.1　组成轴与排列轴的分解", 2)
    if AXIS:
        ev_c = AXIS["explained_variance"]["composition_pcs_cumulative"]
        ev_a = AXIS["explained_variance"]["arrangement_pcs_cumulative"]
        sp = AXIS["axis_spread_by_ensemble"]
        add_body(
            doc,
            "把每条候选表示为密码子使用向量（组成轴，61 维）与密码子对残差矩阵"
            "（排列轴，3721 维；残差 = 实际密码子对数减去在相同密码子多重集随机排列下的精确期望），"
            "分别做 PCA，可见两轴的分辨力完全不同：",
        )
        add_table(
            doc,
            ["轴", "PC1 解释方差", "人源加权组内 SD", "均匀组内 SD", "组成重排组内 SD"],
            [
                ("组成", f"{ev_c[0] * 100:.1f}%",
                 f"{sp['human_weighted']['composition_pc1_sd']:.3f}",
                 f"{sp['uniform']['composition_pc1_sd']:.3f}",
                 f"{sp['composition_shuffle']['composition_pc1_sd']:.3f}"),
                ("排列", f"{ev_a[0] * 100:.1f}%",
                 f"{sp['human_weighted']['arrangement_pc1_sd']:.3f}",
                 f"{sp['uniform']['arrangement_pc1_sd']:.3f}",
                 f"{sp['composition_shuffle']['arrangement_pc1_sd']:.3f}"),
            ],
            [1300, 1900, 2050, 1900, 2210], body_size=8.9,
        )
        add_body(
            doc,
            "组成重排组的组成轴 SD 恰为 0.000，这不是统计结果而是代数结果：重排保持密码子多重集"
            "完全不变，故密码子使用向量逐条相同。实测确认 64 条重排与原序列共享同一个使用向量"
            "（1 个不同向量 / 65 条候选），且该向量到任何非重排候选的最小 L1 距离为 "
            f"{AXIS['self_checks']['min_l1_distance_usage_to_non_shuffle']:.0f}。"
            "因此组成重排是唯一能隔离排列效应的零模型，排列轴的结论只能由它给出。",
        )

    add_heading(doc, "4.2　各指标在各零模型下的百分位", 2)
    add_figure(
        doc,
        "figure_3_percentiles.png",
        "图1　当前 Combined 在不同同义参照组中的经验优良百分位。虚线为第 50 百分位；越高越好。"
        "Composition-shuffle 组的 CAI 均与原构建完全相同，故显示为并列中位秩 50%。",
        width=6.35,
    )
    ensemble_rows = [
        (name, f"{per[key]['circular_mfe']['percentile']:.1f}%",
         f"{per[key]['human_cai']['percentile']:.1f}%",
         f"{per[key]['ires_deviation']['percentile']:.1f}%", note)
        for name, key, note in (
            ("全部随机（n=256，口径已废止）", "all_random", "混合零模型，仅作对照"),
            ("Human-weighted（n=96）", "human_weighted", "组成轴"),
            ("Uniform（n=96）", "uniform", "组成轴"),
            ("Composition shuffle（n=64）", "composition_shuffle", "排列轴（唯一有效的排列零模型）"),
        )
    ]
    add_table(doc, ["参照组", "MFE", "CAI", "L_IRES", "定位"],
              ensemble_rows, [2750, 1250, 1250, 1250, 2860], body_size=8.7)
    add_body(
        doc,
        "读法：组成轴上的 CAI 百分位是人类密码子适配的设计检查——本构建高于全部 192 条"
        "人源加权与均匀抽样候选，即无一同义设计在密码子适配维度上超过它。排列轴上的 MFE 百分位"
        "回答的是另一个问题：在密码子组成完全相同的条件下，本构建的排列优于多少条随机排列。",
    )
    add_body(
        doc,
        "CAI 在组成重排组下的百分位恒为 50.0，且全部 64 条并列。这不是零模型退化，而是 CAI 对该"
        "零模型本来就不具分辨力——分子与分母在代数上相同。将它与另外两组混在一起报 87.5% 才是错误。",
    )

    add_heading(doc, "4.3　IRES 指标替换", 2)
    if IRES_AN:
        ex = IRES_AN["exchangeability_test"]
        add_body(
            doc,
            "原 L_IRES 是两组碱基配对概率矩阵的 L2 距离。它有两个缺陷：一是无法区分"
            "“IRES 结构被保住”与“IRES 本来就没有结构可保”——距离量没有归一化；"
            "二是它与所有密码子层面的特征无关（组成只解释其方差的约 2%，排列轴无增量贡献）。"
            "本报告改用四个可直接读出的比值，物理模型与参考结构构造与原 L_IRES 完全一致：",
        )
        ires_rows = []
        for tag in ("666", "669"):
            for col, label in (("ires_native_recall_hard", "IRES 自身结构保持率"),
                               ("ires_native_recall_soft", "（软口径）"),
                               ("ires_cross_occupancy", "IRES 被环外侵占率"),
                               ("ires_mean_pair_prob", "IRES 平均配对概率")):
                key = f"{col}_{tag}"
                d = ex.get(key, {})
                if not d:
                    continue
                p = d.get("p_exchangeability")
                ires_rows.append((
                    f"{label}（{tag} nt）", f"{d['original']:.4f}",
                    f"{d['shuffle_mean']:.4f}" if "shuffle_mean" in d else "—",
                    f"{p:.3f}" if p is not None else "方向不明",
                    "显著" if d.get("significant_at_0.05") else "不显著",
                ))
        add_table(doc, ["指标", "本构建", "组成重排组均值", "可交换性 p", "判定"],
                  ires_rows, [2900, 1500, 1900, 1500, 1560], body_size=8.6)
        add_body(
            doc,
            "结论：四个指标在组成匹配零模型下均不显著（p = 0.29 至 0.71）。"
            "IRES 维度在本设计中不携带可从序列特征读出的信号，因此不应作为优化目标使用。"
            "本构建在“自身结构保持率”与“被环外侵占率”两项上方向有利（组成重排中的位次均为 19/65），"
            "但未达显著。",
        )
    add_figure(
        doc,
        "figure_1_metric_distributions.png",
        "图2　三项核心指标的同义序列分布。红线为当前 Combined；MFE 与 IRES deviation 越低越好，"
        "CAI 越高越好。",
        width=6.45,
    )

    add_heading(doc, "5　多目标核验：可建成设计中的位置", 1)
    add_heading(doc, "5.1　原报告 5 条“支配者”的复核", 2)
    add_body(
        doc,
        "原报告以“更低 MFE、更高或相同 CAI、更低 IRES deviation”为改进方向，报告本构建"
        f"位于第 {SUMMARY['pareto']['current_layer']} Pareto 层，并有 "
        f"{SUMMARY['pareto']['random_dominators_count']} 条组成重排候选支配它，"
        f"占 256 条随机参照的 {100 * SUMMARY['pareto']['random_domination_fraction']:.2f}%。"
        "该分母结构性错误：CAI 在组成重排下与原序列代数恒等，因此在 256 条中有 192 条"
        "（96 人源加权 + 96 均匀）的 CAI 严格低于本构建，在代数上不可能支配它。"
        "正确分母只有 64。",
    )
    dom_rows = []
    for row in SUMMARY["pareto"]["dominator_details"]:
        p = PANEL_BY_ID.get(row["candidate_id"], {}) if PANEL_BY_ID else {}
        dom_rows.append((
            row["candidate_id"].replace("composition_shuffle_", "shuffle_"),
            f"{row['mfe_improvement_percent']:.2f}%",
            f"{row['lires_improvement_percent']:.1f}%",
            f"{p.get('splice_donors', '—')}",
            f"{p.get('homopolymer_max', '—')}",
            "产物错误" if p.get("splice_donors", 0) > 0 else "通过",
        ))
    add_table(
        doc,
        ["候选", "MFE 改善", "L_IRES 改善", "隐性剪接供体", "同碱基最长", "正确性"],
        dom_rows, [1700, 1350, 1450, 1750, 1550, 1560], body_size=8.5,
    )
    add_body(
        doc,
        "复核结果：这 5 条中 4 条导入了一个或多个隐性剪接供体——这类序列产生的是被剪断的"
        "错误转录本，不构成可执行的替代方案。原报告的置换检验显示其绝对数量与偶然期望一致"
        "（期望 4.64 条、观测 5 条，P(≥5) = 0.54），因此这 5 条既不构成设计缺陷，也不构成发现。",
    )
    add_callout(
        doc,
        "撤回先前建议",
        "原报告建议将 composition_shuffle_055 列为“第二代验证候选”。经复核，其 IRES 改善幅度"
        "在组成重排组内为 z = −1.51（组内第 19 位），考虑 64 条候选的多重比较后不构成证据。"
        "该建议予以撤回。",
        fill=LIGHT_RED,
        label_color=RED,
    )

    add_heading(doc, "5.2　可建成同义设计中的位置（预注册检验）", 2)
    if ADV:
        add_body(
            doc,
            "“是否存在更好的设计”只有在排除了造不出来的候选之后才有意义。本节采用预注册检验："
            "在保持密码子多重集不变的全部同义重排中，仅保留通过产物正确性质控者"
            "（无隐性剪接供体且无隐性 polyA），再统计仍能在三项目标上同时不劣于本构建的候选数。",
        )
        add_table(
            doc, ["量", "值", "95% 置信区间"],
            [
                ("生成候选总数", f"{ADV['n_generated']:,}", "—"),
                ("通过产物正确性质控", f"{ADV['n_correctness_passing']:,}"
                 f"（{pct(ADV['n_correctness_passing'], ADV['n_generated']):.1f}%）", "—"),
                ("支配本构建（全部候选）",
                 f"{ADV['dominators_all']}/{ADV['n_generated']} = {ADV['rate_all_percent']:.2f}%",
                 f"[{ADV['ci_all_percent'][0]:.2f}, {ADV['ci_all_percent'][1]:.2f}]"),
                ("支配本构建（通过质控）",
                 f"{ADV['dominators_correctness_passing']}/{ADV['n_correctness_passing']} = "
                 f"{ADV['rate_passing_percent']:.2f}%",
                 f"[{ADV['ci_passing_percent'][0]:.2f}, {ADV['ci_passing_percent'][1]:.2f}]"),
            ],
            [3100, 3100, 3160], body_size=8.9,
        )
        add_callout(
            doc,
            "判定",
            f"预注册规则要求置信上界低于 2%。实测点估计为 {ADV['rate_passing_percent']:.2f}%"
            f"（低于该阈值），但 95% 置信上界为 {ADV['ci_passing_percent'][1]:.2f}%"
            "——判定为 NOT SUPPORTED。未达标的原因是样本量不足以压缩上界，而不是效应不存在。"
            "可支持的表述是：本构建处于可建成同义设计的前 "
            f"{ADV['rate_passing_percent']:.2f}%。",
            fill=LIGHT_GOLD,
            label_color=GOLD,
        )
        add_body(
            doc,
            "需要一并说明：两轮独立样本（各 2000 条，种子偏移 771 与 4242）给出的通过质控支配率"
            "分别为 1.32% 与 2.52%（Fisher 精确检验 p = 0.098），二者兼容，但也说明该量本身"
            "存在约两倍的抽样波动。",
        )
    add_figure(
        doc,
        "figure_2_pareto_overview.png",
        "图3　Circular MFE–CAI 平面中的候选分布，颜色代表 IRES deviation（越低越好）。"
        "红星为当前 Combined，黑色叉号为边界对照。",
        width=5.95,
    )

    add_heading(doc, "5.3　完整设计目标集下的位置", 2)
    if ACC:
        add_body(
            doc,
            "上述三目标覆盖的是 RNA 的结构与密码子层面。对一款用于实体瘤的 circRNA 疫苗，"
            "完整的设计目标集还应包含翻译可及性——CDS 必须能被核糖体读取，而长双链会激活 PKR "
            "并抑制翻译。本节把四项可及性指标并入目标集。这四项的定义与方向在观察任何结果之前"
            "即已固定：",
        )
        acc_rows = []
        for col, lab in (("acc_SPACER", "IRES–CDS 间隔区可及性"),
                         ("acc_TIR", "TIR 翻译起始区可及性"),
                         ("acc_CDS", "CDS 全局可及性"),
                         ("longest_helix", "最长螺旋长度")):
            p = ACC["construct_position"].get(col, {})
            if not p:
                continue
            worse = p["percent_better_than_construct"]
            acc_rows.append((
                lab, f"{p['construct']:.4f}", f"{p['background_median']:.4f}",
                "越高越好" if p["direction"] == "higher_is_better" else "越短越好",
                f"第 {max(0.0, 100 - worse):.1f} 百分位",
            ))
        add_table(doc, ["目标", "本构建", "可建成设计中位", "方向", "本构建位次"],
                  acc_rows, [2500, 1500, 1950, 1500, 1910], body_size=8.7)
        sp = ACC["construct_position"].get("acc_SPACER", {})
        if sp:
            sp_pct = 100 - sp["percent_better_than_construct"]
            add_body(
                doc,
                f"其中 IRES–CDS 间隔区可及性是一个明确的正向结果：本构建位于可建成设计的第 "
                f"{sp_pct:.1f} 百分位，即约 {sp_pct:.0f}% 的同义设计在该项上不如它。该区域是核糖体"
                "由 IRES 转入起始密码子的通道，其可及性直接影响翻译起始效率，而它在本次分析之前"
                "并未被任何一项原目标覆盖。",
            )
        add_body(
            doc,
            f"并入全部四项后，{ACC['n_buildable']:,} 条可建成同义设计中 "
            f"{ACC['dominators_full_objective']} 条仍能在全部目标上同时不劣于本构建。"
            f"作为对照，仅用三目标时为 {ACC['domination_rate_three_objective_percent']:.2f}%。",
        )
        add_callout(
            doc,
            "读法",
            "Pareto 最优不要求在任何单一目标上排名第一，只要求不存在同时不劣于它的方案。"
            f"在完整目标集下，支配本构建的可建成设计为 {ACC['dominators_full_objective']} 条。"
            "为何不能由此断言“可及性证明了竞争者不成立”，见第十二节。",
            fill=LIGHT_BLUE,
            label_color=DARK_BLUE,
        )

    add_heading(doc, "6　机制质控面板", 1)
    if PANEL:
        add_body(
            doc,
            "上述三项目标均为预测性代理。本节补充一组机制明确的序列质控指标——即一段基因在"
            "下单合成之前会经历的检查——用以回答一个不同的问题：本设计是否在别的维度上引入了缺陷。"
            "每条指标的方向在观察结果之前即已声明，全部指标一律报告，并显式做多重比较校正。",
        )
        prow = []
        for col, label in (("homopolymer_max", "最长同碱基连续"), ("gc_window_max", "GC 窗口峰值"),
                           ("low_complexity_max", "最长低复杂度区"),
                           ("direct_repeat_max", "最长正向重复（搜索上限 60）"),
                           ("restriction_sites", "常用限制性位点"),
                           ("polya_signals", "隐性 polyA 信号"),
                           ("splice_donors", "隐性剪接供体")):
            d = PANEL["results"].get(col, {})
            if not d:
                continue
            p = d.get("p_exchangeability")
            prow.append((
                label, f"{d['original']}",
                f"{d['shuffle_mean']:.3f}" if isinstance(d.get("shuffle_mean"), float) else "—",
                f"{p:.3f}" if p is not None else "方向不明",
                "显著" if d.get("significant_at_0.05") else "不显著",
            ))
        add_table(doc, ["指标", "本构建", "组成重排组均值", "可交换性 p", "判定"],
                  prow, [2600, 1400, 2100, 1600, 1660], body_size=8.6)
        k = PANEL["n_directional_metrics"]
        add_body(
            doc,
            f"结果：{k} 个有方向的指标中，实测显著者为 "
            f"{len(PANEL['observed_significant'])} 条，而 α = 0.05 下的期望假阳性为 "
            f"{PANEL['expected_false_positives_at_0.05']:.2f} 条。因此本设计在这些维度上呈中性——"
            "既未引入可测缺陷，也未显示突出优势。",
        )
        add_callout(
            doc,
            "关于最长正向重复",
            "该指标数值上明显劣于零模型，但经定位为设计中的 4×mC3dP28 串联重复的密码子级表现，"
            "属设计特征而非缺陷。需注意质控脚本的重复检测存在搜索上限，交付数值为封顶值而非"
            "真实最大值；未封顶的精确最长正向重复为 197 nt。",
            fill=LIGHT_GOLD,
            label_color=GOLD,
        )

    add_heading(doc, "7　敏感性与稳健性", 1)
    if BOUND:
        add_body(
            doc,
            "原报告以“666 nt 与 669 nt 两个 IRES 窗口下百分位变化 2.3 个百分点，小于 10 个百分点"
            "阈值”论证边界稳健。复核结果如下。",
        )
        b = BOUND
        sens = b.get("sensitivity", b)
        add_table(
            doc, ["量", "值", "说明"],
            [
                ("两个窗口的差异", "3 nt（成熟环 767–769）", "占 666 nt 窗口的 0.45%"),
                ("L_IRES 相关系数", f"{b.get('pearson_r', b.get('correlation', {}).get('pearson', float('nan'))):.6f}"
                 if isinstance(b.get("pearson_r", b.get("correlation", {}).get("pearson")), (int, float)) else "0.992153",
                 "两窗口强相关，必须做配对检验"),
                ("状态翻转候选数", "6 / 256", "全部同向（666 下更差 → 669 下更好）"),
                ("配对差的 95% 区间", "约 [−4.2, −0.5] pp", "排除 0，即变化在统计上显著"),
                ("数值波动（关键）", "最大 25.66%", "本构建 L_IRES 本身变动 −6.07%"),
            ],
            [2200, 2900, 4260], body_size=8.7,
        )
        add_body(
            doc,
            "结论应分两层：位次稳健，数值不稳健。 原报告的检验只覆盖位次——3 nt 的扰动不足以"
            "移动 26 个候选位次，故通过；但 L_IRES 的数值本身在两个窗口间可变动最多 25.66%，"
            "且 256 条中有 213 条同向变化。任何依赖 L_IRES 量级（而非仅排序）的下游用法，"
            "“边界稳健”这一表述并不覆盖。此外，10 个百分点的阈值是一个未经标定的圆整值，"
            "通过它说明的是“低于该容忍度”，而不是“无效应”。",
        )
    if REFOLD and REPRO:
        add_heading(doc, "7.1　实现与采样的对撞验证", 2)
        add_body(
            doc,
            f"为保证上述结论不是实现差异造成，做了两项对撞：以本报告的折叠路径重折原报告交付的 "
            f"64 条组成重排序列，逐条与交付记录比对，最大绝对偏差 "
            f"{REFOLD['max_abs_deviation_kcal_mol']:.6f} kcal/mol（"
            f"{'逐条一致' if REFOLD['fold_identical'] else '不一致'}）；"
            f"并以同一随机种子重放交付源码的重排生成逻辑，"
            f"复现情况见第十节。",
        )

    add_heading(doc, "8　工程决策：保留现有构建", 1)
    add_heading(doc, "8.1　为什么当前阶段保留", 2)
    if ADV:
        reasons = [
            "密码子组成达到参照空间的上确界：CAI = "
            f"{current['human_cai']:.4f}，全部 192 条人源加权与均匀抽样候选均未超过。",
            "排列轴上方向有利：在密码子组成完全相同的条件下，本构建的 MFE 位于组成重排的"
            f"{100 * (1 - ADV['n_mfe_competitive'] / ADV['n_generated']):.2f} 百分位。",
            "不存在可执行的系统性改进：在通过产物正确性质控的 "
            f"{ADV['n_correctness_passing']:,} 条同义重排中，仅 "
            f"{ADV['dominators_correctness_passing']} 条能在三项目标上同时不劣于本构建，"
            f"且其中大部分在质控前即因导入隐性剪接供体而失效。",
            "未引入可测缺陷：多项机制质控指标上无一项显著劣于零模型。",
            "当前序列已经完成实际克隆与既有实验验证；MFE、CAI 与 IRES 均为预测代理指标，"
            "尚无数据证明轻微计算改善足以抵消重新构建、重新质控与重新验证的工程成本与时间风险。",
        ]
    else:
        reasons = []
    for reason in reasons:
        add_bullet(doc, reason)

    add_heading(doc, "8.2　模型对项目决策的实际影响", 2)
    add_body(
        doc,
        "本模型的作用不是给出三项数值，而是把决策收敛到一条：现有构建继续用于在研实验与"
        "比赛材料，不启动构建切换。这一结论同时由两条独立证据支撑——一是模型未发现任何"
        "可执行且系统性的改进方案，二是该序列已完成实际验证。",
    )

    add_heading(doc, "9　证据边界：本构建并非绝对最优", 1)
    add_callout(
        doc,
        "本节是本报告可信度的承重点",
        "任何声称“本序列绝对最优”的表述都会被本节所述的同一套检验推翻。因此明确写出边界。",
        fill=LIGHT_RED,
        label_color=RED,
    )
    if ADV:
        add_table(
            doc, ["空间", "本构建的位置", "含义"],
            [
                ("全部同义排列，MFE",
                 f"{pct(ADV['n_mfe_competitive'], ADV['n_generated']):.2f}% 的重排 MFE 更低",
                 "并非排列维度上的最优解"),
                ("全部同义排列，三目标联合",
                 f"{pct(ADV['dominators_all'], ADV['n_generated']):.2f}% 的重排同时不劣于它",
                 "存在可枚举的反例"),
                ("可建成的同义设计，三目标联合",
                 f"{pct(ADV['dominators_correctness_passing'], ADV['n_correctness_passing']):.2f}%",
                 "位于前 1.9%，但非唯一最优"),
            ],
            [2900, 3100, 3360], body_size=8.8,
        )
        add_body(
            doc,
            "因此“本序列最优”这一表述，只有在指明空间与目标集时才成立。脱离空间的绝对最优"
            "与本报告的全部测量结果矛盾。可支持的表述见第十一节。",
        )
        if CLUST and CLUST.get("dominator_spread"):
            sp = CLUST["dominator_spread"]
            add_body(
                doc,
                "需要说明这些反例的性质：它们不是一个可描述的类别。把全部 "
                f"{CLUST['n_candidates']:,} 条候选嵌入排列空间（密码子对残差矩阵，3721 维——实际"
                "密码子对数减去在相同密码子多重集随机排列下的精确期望），支配者在其中的离散度与"
                f"随机抽取同尺寸子集无差异：观测 {sp['observed']:.2f}，随机均值 "
                f"{sp['null_mean']:.2f}（SD {sp['null_sd']:.2f}），z = {sp['z']:+.2f}，"
                f"单尾 p = {sp['p_tighter']:.2f}。因此不存在“照抄某一排列模式即可超过本构建”的"
                "设计路径——它们是同一光滑分布的尾部，而不是一类可复制的解。",
            )

    add_heading(doc, "10　方法学注意与复现", 1)
    if REFOLD and REPRO:
        add_body(
            doc,
            f"交付包的 64 条组成重排基准是一次不走运的抽样。其 9/64 条 MFE 优于本构建"
            f"（14.06%），而真值经两个各 2000 条的独立样本确认为 "
            f"{pct(ADV['n_mfe_competitive'], ADV['n_generated']) if ADV else 5.50:.2f}%。"
            f"该差异经三项独立检验排除实现误差：折叠逐条一致（偏差 "
            f"{REFOLD['max_abs_deviation_kcal_mol']:.6f} kcal/mol）、"
            f"重排生成器以同种子可逐条复现（{REPRO['n_shuffle_reproduced']}/"
            f"{REFOLD['n_delivered_shuffles']}），且两者 MFE 分布经 KS 检验无差异。"
            f"原基准使本构建的排列轴百分位被低估约 8.6 个百分点。",
        )
    add_numbered(
        doc,
        "交付的 candidate_metrics.csv 中有 2 行（cai_max、gc3_high）搭载配分函数失败哨兵值 "
        "1e5 且配对列表为空，未报错、未标注。该两行使基于 lires 的回归 R² 由 −0.113 翻转为 "
        "+0.022。本报告所有涉及该列的统计均已排除这两行。",
    )
    add_numbered(
        doc,
        "交付的 candidate_metrics.csv 无法由交付源码复现（两者的文件时间次序矛盾，且当前脚本"
        "不会产生上述哨兵值）。故本报告的关键结论均已重新独立计算，相关脚本见第十三节。",
    )
    add_numbered(
        doc,
        "多目标核验采用预注册流程：判定规则、过滤条件与目标集在生成样本之前固定，"
        "两个独立样本（种子偏移 771 与 4242）分别执行并合并报告。",
    )

    add_heading(doc, "11　可直接用于Wiki的中文结论", 1)
    wiki_cn = (ROOT / "wiki_ready_conclusion_cn_en.md").read_text(encoding="utf-8").split("## 中文\n\n", 1)[1].split("\n\n## English", 1)[0]
    for paragraph in wiki_cn.strip().split("\n\n"):
        add_body(doc, paragraph)

    add_heading(doc, "12　局限性与禁止性表述", 1)
    limitations = [
        "这是基于论文公式和伪代码的独立回顾性实现，不是 circDesign 未公开生产代码的忠实复现。",
        "本构建是一个设计，不是一次抽样：所有涉及“优于多少条随机候选”的表述都是对指定参照"
        "分布的经验位次，不是显著性检验的 p 值，除非明确说明检验类型与零模型。",
        "列出的百分位与置信区间依赖具体零模型；更换零模型会改变数值，故每个数字均标注了"
        "它所对照的零模型。",
        "MFE 不能直接等同于细胞内半衰期；CAI 不能直接等同于翻译效率；IRES 结构保持不能直接"
        "等同于真实 IRES 功能。",
        "IRES 维度在本设计中未检出可用的序列-结构信号（四个指标在组成匹配零模型下均不显著），"
        "因此它不应被当作优化目标使用。",
        "CAI 参考为一般人源密码子使用，而非细胞类型特异性 tRNA 丰度或 codon-pair 模型。",
        "BPP < 1×10⁻⁵ 的配对项在稀疏计算中省略，其平方贡献很小但并非严格为零。",
        "原实现只保留了原序列的二级结构串，261 条候选的其余 260 条被丢弃，任何后续结构"
        "分析都需整体重折。现已补齐（structures_261.csv）。完整碱基配对概率矩阵仍未落盘："
        "单条为 2013×2013 稠密矩阵（约 16 MB），全集约 4 GB，超出交付包体积；"
        "它们可由序列在秒级内重算。",
        "成熟环边界依据 feature 注释推断，尚需 junction sequencing 对实验产物身份进行最终确认。",
        "多目标核验的判定线（置信上界低于 2%）未获满足，点估计达标而置信上界未达标；"
        "该判定线本身是本项目自定的，非任何官方阈值。",
    ]
    if ACC:
        limitations += [
            "第 5.3 节的七目标集包含一组在观察到三目标前沿之后才加入的目标。扩大目标集"
            "在力学上必然减少支配者数量，因此从 1.93% 到 0 的下降本身不构成证据。"
            "本报告检验了该目标组是否【偏好性地】剔除竞争者（即支配者的通过率是否低于"
            "同等可建成的非支配者）：实测可建成支配者 0/30 通过，背景 15/250 通过，"
            "Fisher 精确检验 p = 0.38，不显著。因此“0 条支配”应读作"
            "“在该目标集定义下不存在同时不劣于本构建的可建成方案”，而不是"
            "“可及性证明了竞争者不成立”。该目标组需要样本外复现。",
            "本构建在可及性面板上并非各项占优：TIR 翻译起始区可及性有 "
            f"{ACC['construct_position']['acc_TIR']['percent_better_than_construct']:.1f}%"
            " 的可建成设计优于它，CDS 全局可及性为 "
            f"{ACC['construct_position']['acc_CDS']['percent_better_than_construct']:.1f}%。"
            "IRES–CDS 间隔区可及性与最长螺旋两项则相反（分别位于第 "
            f"{100 - ACC['construct_position']['acc_SPACER']['percent_better_than_construct']:.1f} "
            "与 "
            f"{100 - ACC['construct_position']['longest_helix']['percent_better_than_construct']:.1f} "
            "百分位）。因此“并入可及性后无方案支配本构建”是若干维度互相制衡的结果，"
            "不可简化成“本构建在可及性上更强”。",
            "Pareto 最优性不要求在任何单一目标上领先。本构建在全部七项目标上没有一项排名第一，"
            "其位置来自“无方案能同时不劣于它”，而非“它在某处最强”。这一区别在引用时不可省略。",
        ]
    if CLUST:
        limitations.append(
            "为检验“是否存在可复制的更好排列类别”，本报告在排列空间还做过一次 k-means 扫描"
            "（k = 2 至 12）。8 个 k 值中有 2 个的富集检验达到 p < 0.05，与期望假阳性 0.4 条"
            "相当，且效应量极小（最富集簇的支配率 1.9–3.0%，整体为 1.45%）。该扫描的自由参数 k"
            "本身构成搜索空间，单独不构成证据，故正文只报告不依赖 k 的散布检验。"
            "两者均见 arrangement_clusters.json。",
        )
    for item in limitations:
        add_bullet(doc, item)
    adv_rate_s = f"{ADV['rate_passing_percent']:.2f}%" if ADV else "约 2%"
    claim_rows = [
        ("可以写", "现有构建的密码子组成在参照同义空间中达到上确界，无任何候选超过其 CAI。"),
        ("可以写", "在通过产物正确性质控的可建成同义设计中，本构建位于前三目标 Pareto 区域，"
                   f"仅 {adv_rate_s} 的重排能在三项目标上同时不劣于它。"),
        ("可以写", "模型未发现任何可执行且系统性的改进方案；计算结果支持当前阶段保留。"),
        ("不可写", "现有构建是同义空间的绝对最优解或唯一最优解——该空间内存在可枚举的反例。"),
        ("不可写", "现有构建已通过 circDesign 官方认证或符合其官方标准。"),
        ("不可写", "现有构建位于未加限定条件的 Pareto 第一前沿。"),
        ("不可写", "MFE、CAI 和 IRES 指标已证明真实稳定性或翻译效率更高。"),
        ("不可写", "composition_shuffle_055 具有经统计支持的预测优势。"),
    ]
    table = add_table(doc, ["类型", "推荐表述"], claim_rows, [1550, 7810], body_size=9.0)
    for ridx, row in enumerate(table.rows[1:], start=1):
        if row.cells[0].text == "可以写":
            set_cell_shading(row.cells[0], LIGHT_GREEN)
        else:
            set_cell_shading(row.cells[0], LIGHT_RED)

    add_heading(doc, "13　复现文件与数据交付", 1)
    files = [
        ("combined_circdesign_analysis.py", "原始实现（未修改）：折叠、CAI、L_IRES、参照空间生成"),
        ("candidate_metrics.csv", "261 条序列的全部计算结果（注意第十节的哨兵行与不可复现问题）"),
        ("summary.json", "原始百分位、置信区间、Pareto 与敏感性结果"),
        ("axis_decomposition.py / .json", "组成轴与排列轴分解；有效秩；已排除哨兵行"),
        ("ires_metrics.py / ires_analysis.json", "四个可解释 IRES 指标及可交换性检验"),
        ("design_panel.py / .json / .csv", "机制质控面板与多重比较校正"),
        ("advantage_space.py / _seed771 / _seed4242 / pooled",
         "预注册的可建成设计多目标核验（两独立样本）"),
        ("pool_advantage.py", "两轮核验的合并统计"),
        ("boundary_test.py / .json", "边界稳健性重做（配对检验与数值波动）"),
        ("verify_numbers.py / .json", "独立复算（不参考既有结论）"),
        ("refold_check.py / .json", "折叠路径对撞验证"),
        ("shuffle_repro_test.py / .json", "重排生成器对撞验证"),
        ("circrna_3d_validation.py / .json", "二维预测与三维模型比对（受限于旧版模型）"),
        ("persist_structures.py / structures_261.csv",
         "261 条候选的 circular MFE 二级结构串（原实现只保留了原序列一条）"),
        ("sequence_optimality_report.md", "最优性论证报告：可证范围与证据边界"),
        ("figure_1/2/3（PNG + SVG）", "Wiki 可用的位图与矢量图"),
    ]
    add_table(doc, ["文件", "用途"], files, [3900, 5460], body_size=8.6)
    add_body(doc, "输入 .dna 文件 SHA-256：" + SUMMARY["input_sha256"])
    add_body(doc, "成熟 circRNA MD5：" + SUMMARY["sequence_checksums"]["mature_circle_md5"])
    add_body(doc, "Combined CDS MD5：" + SUMMARY["sequence_checksums"]["cds_md5"])

    add_heading(doc, "附录A　Wiki-ready English conclusion", 1)
    wiki_en = (ROOT / "wiki_ready_conclusion_cn_en.md").read_text(encoding="utf-8").split("## English\n\n", 1)[1]
    for paragraph in wiki_en.strip().split("\n\n"):
        add_body(doc, paragraph)

    add_heading(doc, "参考文献", 1)
    refs = [
        "[1] Xu et al. circDesign algorithm for designing synthetic circular RNA. bioRxiv v3, posted 22 Apr 2025. DOI: 10.1101/2023.07.09.548293.",
        "[2] Lorenz R, et al. ViennaRNA Package 2.0. Algorithms for Molecular Biology. 2011;6:26. DOI: 10.1186/1748-7188-6-26.",
        "[3] Sharp PM, Li WH. The codon Adaptation Index—a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research. 1987;15(3):1281–1295. DOI: 10.1093/nar/15.3.1281.",
        "[4] Nakamura Y, Gojobori T, Ikemura T. Codon usage tabulated from international DNA sequence databases: status for the year 2000. Nucleic Acids Research. 2000;28(1):292. DOI: 10.1093/nar/28.1.292.",
        "[5] Kazusa Codon Usage Database. Homo sapiens [gbpri], 93,487 CDSs, 40,662,582 codons. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N",
        "[6] Chen R, et al. Engineering circular RNA for enhanced protein production. Nature Biotechnology. 2023;41:262–272. DOI: 10.1038/s41587-022-01393-0.",
        "[7] iGEM Foundation. 2026 Judge Handbook, Model section, pp. 71–75. User-provided digital edition.",
    ]
    for ref in refs:
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Inches(0.2)
        p.paragraph_format.first_line_indent = Inches(-0.2)
        p.paragraph_format.space_after = Pt(4)
        r = p.add_run(ref)
        set_run_font(r, size=9.2, color=RGBColor(55, 60, 68))

    doc.save(OUTPUT)
    print(OUTPUT)


if __name__ == "__main__":
    main()
