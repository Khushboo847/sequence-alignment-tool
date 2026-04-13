from flask import Flask, render_template_string, request

app = Flask(__name__)

SIMILAR_GROUPS = [
    set("AVLIM"),
    set("KRH"),
    set("DE"),
    set("STNQ"),
    set("FYW"),
]

def are_similar(a, b):
    return any(a in g and b in g for g in SIMILAR_GROUPS)


def classify(a, b):
    if a == b:
        return "|"
    elif a != "-" and b != "-" and are_similar(a, b):
        return ":"
    elif a != "-" and b != "-":
        return "."
    return " "


def align(seq1, seq2):
    num1, num2 = 0, 0

    pos1_line = []
    seq1_line = []
    mid_line = []
    seq2_line = []
    pos2_line = []

    for a, b in zip(seq1, seq2):
        if a != "-":
            num1 += 1
            pos1_line.append(f"{num1:^6}")
        else:
            pos1_line.append(" " * 6)
        seq1_line.append(f"{a:^6}")

        if b != "-":
            num2 += 1
            pos2_line.append(f"{num2:^6}")
        else:
            pos2_line.append(" " * 6)
        seq2_line.append(f"{b:^6}")

        mid_line.append(f"{classify(a, b):^6}")

    return {
        "pos1": "".join(pos1_line),
        "seq1": "".join(seq1_line),
        "mid": "".join(mid_line),
        "seq2": "".join(seq2_line),
        "pos2": "".join(pos2_line),
    }


HTML = """
<!DOCTYPE html>
<html>
<head>
    <title>Sequence Alignment Tool</title>
    <style>
        body { font-family: monospace; padding: 20px; }
        input, button { padding: 10px; margin: 5px; width: 300px; }
        .output { margin-top: 20px; white-space: pre; }
    </style>
</head>
<body>
    <h2>Sequence Alignment Tool</h2>
    <form method="post">
        <input type="text" name="seq1" placeholder="Sequence 1" required><br>
        <input type="text" name="seq2" placeholder="Sequence 2" required><br>
        <button type="submit">Align</button>
    </form>

    {% if result %}
    <div class="output">
Pos1: {{ result.pos1 }}
Seq1: {{ result.seq1 }}
      {{ result.mid }}
Seq2: {{ result.seq2 }}
Pos2: {{ result.pos2 }}
    </div>
    {% endif %}
</body>
</html>
"""

@app.route("/", methods=["GET", "POST"])
def home():
    result = None
    if request.method == "POST":
        seq1 = request.form["seq1"].upper()
        seq2 = request.form["seq2"].upper()
        result = align(seq1, seq2)
    return render_template_string(HTML, result=result)


if __name__ == "__main__":
    app.run()
