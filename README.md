#  Sequence Alignment Visualizer

A lightweight and user-friendly web tool for visualizing sequence alignments with clear residue comparison and automatic numbering.

---

##  Features

-  **Residue-wise comparison**
  - Identical (`|`)
  - Similar (`:`)
  - Mismatch (`.`)

-  **Automatic residue numbering**
  - Independent numbering for both sequences
  - Gaps (`-`) are skipped correctly

-  **Improved readability**
  - Clean, aligned output
  - No need for manual counting

-  **Lightweight & fast**
  - No heavy bioinformatics libraries required

---

##  How to run

link : https://sequence-alignment-tool.onrender.com

---

##  Tech Stack

- Python
- Flask
- HTML/CSS

---

##  Example Output

Pos1:    1      2      3             4      5
Seq1:    A      L      V      -      K      T
         |      :      |      .      |
Seq2:    A      I      V      R      -      T
Pos2:    1      2      3      4             5

---


##  Motivation

Traditional sequence alignment tools often make it difficult to track residue positions, especially in the presence of gaps.

This tool focuses on:
- Eliminating manual counting
- Improving clarity
- Making alignment interpretation faster and easier

---

##  Future Improvements

-  Color-coded alignment (like BLAST)
-  Identity & similarity statistics
-  Export results (TXT/CSV)
-  Support for longer sequences with scrolling view

---

##  Contributing

Contributions are welcome! Feel free to fork the repo and submit a pull request.

---

##  License

This project is open-source and available under the MIT License.

---

##  Author

Developed as a bioinformatics utility tool to simplify sequence alignment visualization.
