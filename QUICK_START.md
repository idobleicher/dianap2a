# Quick Start Guide - Viral 2A Protein Finder

## What You Have

A complete Python tool for identifying viral 2A proteins from ribosomal skipping!

## Files in This Project

- **viral_2a_finder.py** - Main analysis script
- **README.md** - Detailed documentation
- **QUICK_START.md** - This file (quick reference)

## How to Run

```bash
python viral_2a_finder.py
```

or

```bash
py viral_2a_finder.py
```

That's it! The script will:

1. Download 500 viral polyproteins from UniProt
2. Find 2A motifs (DxExNPGP patterns)
3. Extract downstream proteins starting with Pro (P)
4. Classify by viral family
5. Generate output files automatically

## Output Files

After running, you'll get 4 files (with timestamps):

1. **2a_proteins_summary_TIMESTAMP.csv**

   - Quick summary table
   - Open in Excel/Google Sheets

2. **2a_proteins_detailed_TIMESTAMP.xlsx**

   - Full Excel workbook
   - Multiple sheets: All Results, High Confidence, By Family

3. **2a_downstream_proteins_TIMESTAMP.fasta**

   - All downstream protein sequences
   - Ready for BLAST or further analysis

4. **2a_downstream_HIGH_CONFIDENCE_TIMESTAMP.fasta**
   - Only validated 2A proteins from known families
   - Best for publication-quality results

## Current Results (from your first run)

**Found: 16 viral 2A proteins (all high confidence!)**

### By Family

- **Picornaviridae**: 13 proteins
  - 10x Foot-and-mouth disease virus (various serotypes)
  - Seneca Valley virus
  - Ljunganvirus
- **Picorna-like viruses**: 3 proteins
  - Acute bee paralysis virus
  - Ectropis obliqua picorna-like virus (2 motifs)

### 2A Motif Patterns Found

- **DVESNPGP** - 12 occurrences (most common)
- **DVETNPGP** - 2 occurrences
- **DIETNPGP** - 1 occurrence
- **DIESNPGP** - 1 occurrence

## Customization

### Get More Results

Edit `viral_2a_finder.py`, line 324:

```python
finder.run_analysis(batch_size=1000)  # Change from 500 to 1000+
```

### Add More Virus Families

Edit the `KNOWN_2A_FAMILIES` dictionary (lines 25-50) to add your own virus patterns.

## Key Features

- [x] Finds canonical 2A motifs: **D[VIL]E[ST]NPGP**
- [x] Validates downstream proteins start with **Proline (P)**
- [x] Ribosomal skipping at **NPG|P junction**
- [x] Excludes protease cleavage sites
- [x] High-confidence classification
- [x] Publication-ready FASTA files

## Scientific Notes

### What is 2A?

2A is a short peptide sequence (~18-20 amino acids) that causes **ribosomal skipping**:

- During translation, ribosome reaches NPG|P junction
- Peptide bond formation is **skipped**
- Two separate proteins are produced
- Downstream protein begins with **Proline**

### Why It Matters

Used by:

- Picornaviruses (FMDV, poliovirus, etc.)
- Biotechnology (multi-gene expression)
- Gene therapy vectors
- Synthetic biology constructs

### Classic Example: FMDV

Foot-and-mouth disease virus (FMDV) 2A:

- Sequence: **LLNFDLLKLAGDVESNPGP**
- Cleavage: ...NPG|P...
- ~18-19 residues long
- Most studied and widely used

## Next Steps

### Analyze Your Results

1. Open the Excel file to explore data
2. Check the FASTA files for sequences
3. Run BLAST on interesting hits
4. Compare motif patterns

### Run More Comprehensive Search

```bash
# Edit viral_2a_finder.py
# Change batch_size to 2000 or 5000
py -3.13 viral_2a_finder.py
```

### Focus on Specific Families

You can modify the UniProt query in the script to target specific virus families.

## Common Questions

**Q: Why only 16 results from 500 sequences?**
A: 2A motifs are relatively rare. Most viral polyproteins use proteases instead.

**Q: Can I search for different motifs?**
A: Yes! Edit the regex patterns in the `find_2a_motifs()` function.

**Q: Are these true 2A proteins?**
A: All 16 hits are from known 2A-containing families, making them high-confidence results!

**Q: What about other cleavage mechanisms?**
A: This tool specifically finds 2A ribosomal skipping, not protease cleavage.

## Citation

If you use this tool in research, please cite:

- UniProt database
- Biopython
- Relevant 2A literature (see README.md)

## Troubleshooting

**Script fails to download:**

- Check internet connection
- UniProt API might be busy
- Try again in a few minutes

**Want more hits:**

- Increase batch_size
- Try broader regex patterns
- Search specific viral families

## Resources

- UniProt: https://www.uniprot.org
- Biopython: https://biopython.org
- 2A peptides review: Ryan & Drew (1994) EMBO J

---

**You're all set! Your Python 3.13.11 installation is working perfectly.**

Happy viral protein hunting! 🦠🔬
