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
    from itertools import zip_longest

    num1, num2 = 0, 0
    blocks = []

    for a, b in zip_longest(seq1, seq2, fillvalue='-'):
        if a != "-":
            num1 += 1
            p1 = str(num1)
        else:
            p1 = ""

        if b != "-":
            num2 += 1
            p2 = str(num2)
        else:
            p2 = ""

        match = classify(a, b)

        blocks.append({
            "a": a,
            "b": b,
            "p1": p1,
            "p2": p2,
            "match": match
        })

    return blocks


HTML = """
<!DOCTYPE html>
<html>
<head>
    <title>Sequence Alignment Visualizer</title>
    <style>
        body {
            font-family: Arial;
            background: #0f172a;
            color: white;
            text-align: center;
            padding: 20px;
        }

        h1 {
            margin-bottom: 20px;
        }

        input {
            padding: 12px;
            margin: 10px;
            width: 300px;
            border-radius: 8px;
            border: none;
        }

        button {
            padding: 12px 20px;
            border: none;
            border-radius: 8px;
            background: #22c55e;
            color: black;
            font-weight: bold;
            cursor: pointer;
        }

        .alignment {
            margin-top: 30px;
            display: inline-block;
        }

        .row {
            display: flex;
            justify-content: center;
        }

        .cell {
            width: 50px;
            height: 50px;
            margin: 2px;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            border-radius: 6px;
            font-weight: bold;
        }

        .identical { background: #22c55e; }
        .similar { background: #eab308; color: black; }
        .mismatch { background: #ef4444; }
        .gap { background: #374151; }

        .pos {
            font-size: 10px;
            opacity: 0.7;
        }

        .res {
            font-size: 16px;
        }
    </style>
</head>
<body>

<h1>🧬 Sequence Alignment Visualizer</h1>

<form method="post">
    <input type="text" name="seq1" placeholder="Enter Sequence 1" required><br>
    <input type="text" name="seq2" placeholder="Enter Sequence 2" required><br>
    <button type="submit">Align</button>
</form>

{% if blocks %}
<div class="alignment">

    <!-- Seq1 -->
    <div class="row">
        {% for b in blocks %}
        <div class="cell">
            <div class="pos">{{ b.p1 }}</div>
            <div class="res">{{ b.a }}</div>
        </div>
        {% endfor %}
    </div>

    <!-- Match -->
    <div class="row">
        {% for b in blocks %}
        {% set cls = 'gap' %}
        {% if b.match == '|' %}
            {% set cls = 'identical' %}
        {% elif b.match == ':' %}
            {% set cls = 'similar' %}
        {% elif b.match == '.' %}
            {% set cls = 'mismatch' %}
        {% endif %}
        <div class="cell {{ cls }}">
            {{ b.match }}
        </div>
        {% endfor %}
    </div>

    <!-- Seq2 -->
    <div class="row">
        {% for b in blocks %}
        <div class="cell">
            <div class="pos">{{ b.p2 }}</div>
            <div class="res">{{ b.b }}</div>
        </div>
        {% endfor %}
    </div>

</div>
{% endif %}

</body>
</html>
"""

@app.route("/", methods=["GET", "POST"])
def home():
    blocks = None
    if request.method == "POST":
        seq1 = request.form["seq1"].upper()
        seq2 = request.form["seq2"].upper()
        blocks = align(seq1, seq2)
    return render_template_string(HTML, blocks=blocks)


if __name__ == "__main__":
    app.run()
