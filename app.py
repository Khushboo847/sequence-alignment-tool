from flask import Flask, render_template_string, request, send_file
import io

app = Flask(__name__)
VISIBLE_CLASSES = ("identical", "strongly-similar", "weakly-similar")

# Conservative substitutions.
SIMILAR_GROUPS = [
    set("AVLIM"),
    set("KRH"),
    set("DE"),
    set("STNQ"),
    set("FYW"),
]

# Broader families used as semi-similar substitutions.
SEMI_SIMILAR_GROUPS = [
    set("AG"),
    set("ILMV"),
    set("FYW"),
    set("STNQ"),
    set("KRH"),
    set("DE"),
    set("CP"),
]


def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2):
    """
    Needleman-Wunsch algorithm for global pairwise alignment.
    Returns aligned seq1 and seq2.
    """
    n = len(seq1)
    m = len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap_penalty
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (
                match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            )
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty
            score[i][j] = max(diag, up, left)

    align1 = []
    align2 = []
    i, j = n, m

    while i > 0 or j > 0:
        if (
            i > 0
            and j > 0
            and score[i][j]
            == score[i - 1][j - 1]
            + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
        ):
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i - 1][j] + gap_penalty:
            align1.append(seq1[i - 1])
            align2.append("-")
            i -= 1
        else:
            align1.append("-")
            align2.append(seq2[j - 1])
            j -= 1

    align1.reverse()
    align2.reverse()
    return "".join(align1), "".join(align2)


def in_same_group(a, b, groups):
    return any(a in group and b in group for group in groups)


def classify_residue(a, b):
    if a == "-" or b == "-":
        return "other", " "
    if a == b:
        return "identical", "|"
    if in_same_group(a, b, SIMILAR_GROUPS):
        return "strongly-similar", ":"
    if in_same_group(a, b, SEMI_SIMILAR_GROUPS):
        return "weakly-similar", ";"
    return "other", " "


def build_alignment_report(aligned1, aligned2):
    counts = {
        "identical": 0,
        "strongly-similar": 0,
        "weakly-similar": 0,
    }

    rows = []
    symbol_line = []
    seq1_pos_line = []
    seq2_pos_line = []
    seq1_pos = 0
    seq2_pos = 0

    for aln_pos, (a, b) in enumerate(zip(aligned1, aligned2), start=1):
        pos1 = "-"
        pos2 = "-"

        if a != "-":
            seq1_pos += 1
            pos1 = seq1_pos
        if b != "-":
            seq2_pos += 1
            pos2 = seq2_pos

        category, symbol = classify_residue(a, b)
        if category in counts:
            counts[category] += 1
        symbol_line.append(symbol)
        seq1_pos_line.append(str(pos1))
        seq2_pos_line.append(str(pos2))
        rows.append(
            {
                "aln_pos": aln_pos,
                "seq1_pos": pos1,
                "seq1_res": a,
                "seq2_pos": pos2,
                "seq2_res": b,
                "category": category,
            }
        )

    visible_rows = [row for row in rows if row["category"] in VISIBLE_CLASSES]

    return {
        "aligned_seq1": aligned1,
        "aligned_seq2": aligned2,
        "symbol_line": "".join(symbol_line),
        "seq1_pos_line": seq1_pos_line,
        "seq2_pos_line": seq2_pos_line,
        "counts": counts,
        "rows": visible_rows,
    }


def count_non_gap(text):
    return sum(1 for ch in text if ch != "-")


def normalize_sequence(raw):
    """
    Normalize user input sequence:
    - supports pasted FASTA/multiline input
    - removes whitespace, numbers, and non-residue symbols
    - keeps letters and gaps only
    """
    residues = []
    for line in raw.splitlines():
        line = line.strip().upper()
        if not line or line.startswith(">"):
            continue
        for ch in line:
            if ch.isalpha() or ch == "-":
                residues.append(ch)
    return "".join(residues)


def format_clustal_like_blocks(seq1_name, seq2_name, report, block_size=60):
    aligned1 = report["aligned_seq1"]
    aligned2 = report["aligned_seq2"]
    marks = report["symbol_line"]
    lines = []
    pos1 = 0
    pos2 = 0
    name_width = max(len(seq1_name), len(seq2_name), 6) + 2

    for start in range(0, len(aligned1), block_size):
        part1 = aligned1[start : start + block_size]
        part2 = aligned2[start : start + block_size]
        partm = marks[start : start + block_size]

        pos1 += count_non_gap(part1)
        pos2 += count_non_gap(part2)

        lines.append(f"{seq1_name:<{name_width}}{part1}  {pos1}")
        lines.append(f"{'':<{name_width}}{partm}")
        lines.append(f"{seq2_name:<{name_width}}{part2}  {pos2}")
        lines.append("")

    return lines


HTML = """
<!DOCTYPE html>
<html>
<head>
<title>Sequence Alignment Analyzer</title>
<style>
body {
    font-family: monospace;
    background: #0f172a;
    color: white;
    text-align: center;
    padding: 30px;
}

.sequence-box {
    padding: 10px;
    margin: 10px;
    width: min(760px, 90vw);
    height: 60px;
    border-radius: 6px;
    border: none;
}

button {
    padding: 10px 20px;
    border-radius: 6px;
    border: none;
    background: #22c55e;
    cursor: pointer;
    font-weight: 700;
}

.summary {
    margin-top: 24px;
}

.top-right-download {
    position: fixed;
    top: 16px;
    right: 16px;
    z-index: 10;
}

.summary table,
.alignment,
.report {
    width: min(1100px, 95vw);
}

.summary table {
    margin: 0 auto;
    border-collapse: collapse;
}

.summary td, .summary th {
    border: 1px solid #333;
    padding: 8px 10px;
}

.alignment {
    margin-top: 30px;
    font-size: 15px;
    line-height: 1.4;
    border-collapse: collapse;
    margin-left: auto;
    margin-right: auto;
}

.alignment td {
    padding: 5px;
    text-align: center;
    border: 1px solid #333;
}

.pos-row td {
    color: #93c5fd;
    font-size: 12px;
}

.report {
    margin: 18px auto 0 auto;
    border-collapse: collapse;
}

.report th, .report td {
    border: 1px solid #333;
    padding: 6px 8px;
    text-align: center;
}

.identical {
    background-color: #4caf50;
    color: white;
}

.strongly-similar {
    background-color: #ffeb3b;
    color: black;
}

.weakly-similar {
    background-color: #ff9800;
    color: black;
}

.download {
    margin-top: 20px;
}

.muted {
    color: #cbd5e1;
}

label {
    display: block;
    font-weight: 700;
    margin-top: 10px;
}
</style>
</head>
<body>

<h1>Sequence Alignment Analyzer</h1>
<p class="muted">Enter two sequences in separate boxes to get counts and residue mapping.</p>

<form method="post">
    <label for="seq1">Sequence 1</label>
    <textarea id="seq1" class="sequence-box" name="seq1" placeholder="Enter Sequence 1" required>{{ seq1_text }}</textarea><br>
    <label for="seq2">Sequence 2</label>
    <textarea id="seq2" class="sequence-box" name="seq2" placeholder="Enter Sequence 2" required>{{ seq2_text }}</textarea><br>
    <button type="submit">Align</button>
</form>

{% if error_message %}
<p style="color:#fca5a5">{{ error_message }}</p>
{% endif %}

{% if report %}
<form method="post" action="/download" class="top-right-download">
    <input type="hidden" name="seq1" value="{{ seq1_text }}">
    <input type="hidden" name="seq2" value="{{ seq2_text }}">
    <button type="submit">Download Output</button>
</form>

<div class="summary">
    <h2>Counts</h2>
    <table>
        <tr>
            <th>Identical</th>
            <th>Strongly Similar</th>
            <th>Weakly Similar</th>
        </tr>
        <tr>
            <td class="identical">{{ report.counts["identical"] }}</td>
            <td class="strongly-similar">{{ report.counts["strongly-similar"] }}</td>
            <td class="weakly-similar">{{ report.counts["weakly-similar"] }}</td>
        </tr>
    </table>
</div>

<table class="alignment">
    <tr>
        <td>Aln Pos</td>
        {% for i in range(1, report.aligned_seq1|length + 1) %}
        <td>{{ i }}</td>
        {% endfor %}
    </tr>
    <tr class="pos-row">
        <td>Seq1 Pos</td>
        {% for p in report.seq1_pos_line %}
        <td>{{ p }}</td>
        {% endfor %}
    </tr>
    <tr>
        <td>Seq1</td>
        {% for char in report.aligned_seq1 %}
        <td>{{ char }}</td>
        {% endfor %}
    </tr>
    <tr>
        <td>Mark</td>
        {% for char in report.symbol_line %}
        <td>{{ char }}</td>
        {% endfor %}
    </tr>
    <tr>
        <td>Seq2</td>
        {% for char in report.aligned_seq2 %}
        <td>{{ char }}</td>
        {% endfor %}
    </tr>
    <tr class="pos-row">
        <td>Seq2 Pos</td>
        {% for p in report.seq2_pos_line %}
        <td>{{ p }}</td>
        {% endfor %}
    </tr>
</table>

<h2>Residue Mapping</h2>
<table class="report">
    <tr>
        <th>Aln Pos</th>
        <th>Seq1 Pos</th>
        <th>Seq1 Res</th>
        <th>Seq2 Pos</th>
        <th>Seq2 Res</th>
        <th>Class</th>
    </tr>
    {% for row in report.rows %}
    <tr class="{{ row.category }}">
        <td>{{ row.aln_pos }}</td>
        <td>{{ row.seq1_pos }}</td>
        <td>{{ row.seq1_res }}</td>
        <td>{{ row.seq2_pos }}</td>
        <td>{{ row.seq2_res }}</td>
        <td>{{ row.category }}</td>
    </tr>
    {% endfor %}
</table>
{% endif %}

</body>
</html>
"""


@app.route("/", methods=["GET", "POST"])
def home():
    report = None
    error_message = ""
    seq1_text = request.form.get("seq1", "")
    seq2_text = request.form.get("seq2", "")

    if request.method == "POST":
        seq1_text = request.form.get("seq1", "")
        seq2_text = request.form.get("seq2", "")
        seq1 = normalize_sequence(seq1_text)
        seq2 = normalize_sequence(seq2_text)

        if not seq1 or not seq2:
            error_message = "Please provide both Sequence 1 and Sequence 2."
        else:
            aligned1, aligned2 = needleman_wunsch(seq1, seq2)
            report = build_alignment_report(aligned1, aligned2)

    return render_template_string(
        HTML,
        report=report,
        error_message=error_message,
        seq1_text=seq1_text,
        seq2_text=seq2_text,
    )


@app.route("/download", methods=["POST"])
def download():
    seq1 = normalize_sequence(request.form.get("seq1", ""))
    seq2 = normalize_sequence(request.form.get("seq2", ""))

    if seq1 and seq2:
        aligned1, aligned2 = needleman_wunsch(seq1, seq2)
        report = build_alignment_report(aligned1, aligned2)

        output_lines = [
            "CLUSTAL-LIKE PAIRWISE ALIGNMENT",
            "",
            "Legend: | identical, : strongly similar, ; weakly similar",
            "",
        ]

        output_lines.extend(format_clustal_like_blocks("Seq1", "Seq2", report))
        output_lines.extend(
            [
            "Counts",
            f"Identical: {report['counts']['identical']}",
            f"Strongly Similar: {report['counts']['strongly-similar']}",
            f"Weakly Similar: {report['counts']['weakly-similar']}",
            "",
            "Residue Mapping (identical/strongly-similar/weakly-similar only)",
            "AlnPos\tSeq1Pos\tSeq1Res\tSeq2Pos\tSeq2Res\tClass",
            ]
        )

        for row in report["rows"]:
            output_lines.append(
                f"{row['aln_pos']}\t{row['seq1_pos']}\t{row['seq1_res']}\t"
                f"{row['seq2_pos']}\t{row['seq2_res']}\t{row['category']}"
            )

        output = "\n".join(output_lines)
    else:
        output = "Please provide both Sequence 1 and Sequence 2."

    buffer = io.BytesIO()
    buffer.write(output.encode())
    buffer.seek(0)

    return send_file(
        buffer,
        as_attachment=True,
        download_name="alignment.txt",
        mimetype="text/plain",
    )


if __name__ == "__main__":
    app.run(debug=True)
