# TASSEL Cophenetic Correlation Calculator

A small Streamlit app that tests how well a dendrogram drawn in TASSEL
represents the genetic distance matrix it was built from.

TASSEL calculates the distance matrix and the tree, but it does not report
the cophenetic correlation or a Mantel test. This app takes the two TASSEL
output files and calculates them.

## Input

1. **Distance matrix** exported from TASSEL
   (*Analysis > Relatedness > Distance Matrix*, then *File > Save As*).
2. **Tree in Newick format** exported from TASSEL
   (*Analysis > Relatedness > Cladogram*, then *File > Save As*).

Both files must contain the same taxa. Branch lengths written with a decimal
comma (as TASSEL does on computers with a comma decimal separator) are
accepted.

## What it calculates

- **Cophenetic matrix**: the distance between every pair of taxa on the tree
  (the sum of the branch lengths on the path that connects them).
- **Cophenetic correlation (r)**: Pearson correlation between the distance
  matrix and the cophenetic matrix.
- **Mantel test**: the taxa of the cophenetic matrix are randomly permuted
  (default 9,999 times); *P* is the proportion of permutations with a
  correlation equal to or higher than the observed one.

## Output

- r and *P* on the screen
- the tree as drawn from the uploaded Newick file
- the cophenetic matrix (CSV)
- a results table (CSV)

## How to run

```
git clone https://github.com/ibrahimcelikmbg15/TASSEL-Cophenetic-Calculator.git
cd TASSEL-Cophenetic-Calculator
pip install -r requirements.txt
streamlit run tassel_cophenetic_calculator.py
```

## Example

For 35 carob (*Ceratonia siliqua* L.) genotypes and 32,848 SNP loci, the UPGMA
dendrogram had a cophenetic correlation of r = 0.8965 (*P* = 0.0001, 9,999
permutations, 595 pairs).
