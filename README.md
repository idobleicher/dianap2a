# Viral 2A Protein Finder

A Python tool to identify viral proteins containing functional 2A motifs that undergo ribosomal skipping.

## What are 2A proteins?

2A proteins are short viral peptide sequences that cause ribosomal skipping during translation. This mechanism allows viruses to produce multiple proteins from a single open reading frame (ORF) without true proteolytic cleavage. The hallmark of 2A-mediated skipping is the characteristic motif **DxExNPGP**, where the ribosome "skips" the formation of a peptide bond between the glycine (G) and proline (P).

## Features

- [x] Downloads viral polyproteins from UniProt
- [x] Scans for canonical 2A motifs (D[VIL]E[ST]NPGP)
- [x] Extracts downstream protein sequences starting with Proline
- [x] Filters by known 2A-containing viral families
- [x] Generates comprehensive output files (CSV, Excel, FASTA)
- [x] Provides high-confidence vs low-confidence classifications

## Requirements

- Python 3.13+
- Biopython
- Requests
- Pandas
- OpenPyXL

## Installation

```bash
pip install biopython requests pandas openpyxl
```

## Usage

### Basic Usage

```bash
python viral_2a_finder.py
```

or

```bash
py viral_2a_finder.py
```

This will:
1. Download up to 500 viral polyprotein sequences from UniProt
2. Scan for 2A motifs
3. Extract downstream proteins
4. Generate output files

### Customizing the Search

Edit the `main()` function to change the batch size:

```python
finder.run_analysis(batch_size=1000)  # Download more sequences
```

## Output Files

The script generates multiple output files with timestamps:

1. **`2a_proteins_summary_TIMESTAMP.csv`**
   - Quick overview of all results
   - Columns: virus name, accession, family, 2A motif, position, downstream length, confidence

2. **`2a_proteins_detailed_TIMESTAMP.xlsx`**
   - Comprehensive Excel workbook with multiple sheets:
     - All Results
     - High Confidence (known 2A families)
     - By Family (summary statistics)

3. **`2a_downstream_proteins_TIMESTAMP.fasta`**
   - FASTA file with all downstream protein sequences
   - Sequences start with Proline (P)

4. **`2a_downstream_HIGH_CONFIDENCE_TIMESTAMP.fasta`**
   - Only high-confidence hits from known 2A families

## Known 2A-Containing Viral Families

The script prioritizes these families (high confidence):

- **Picornaviridae** (foot-and-mouth disease, poliovirus)
- **Dicistroviridae** (cricket paralysis virus)
- **Iflaviridae** (insect viruses)
- Teschovirus
- Aphthovirus
- Enterovirus
- Rhinovirus
- Cardiovirus
- Kobuvirus

## 2A Motif Patterns

The script searches for two patterns:

1. **Canonical**: `D[VIL]E[ST]NPGP`
   - Most common and well-validated
   - Examples: DVESNPGP, DIESTNPGP

2. **Broader**: `D..E.NPGP`
   - More permissive pattern
   - Captures less common variants

## Understanding the Results

### Confidence Levels

- **HIGH**: Virus belongs to a known 2A-containing family
- **LOW**: 2A motif found but not in known family (may be false positive)

### Key Columns

- **2a_motif**: The exact sequence of the 2A peptide
- **2a_position**: Position in the polyprotein
- **downstream_length**: Length of the downstream protein
- **downstream_sequence**: Full sequence starting with P

## Scientific Background

### Ribosomal Skipping Mechanism

1. Translation proceeds normally until the 2A motif
2. At the NPG↓P junction, the ribosome skips peptide bond formation
3. The downstream protein is released with Proline at the N-terminus
4. Translation continues to produce a separate protein

### vs Protease Cleavage

This tool focuses on **ribosomal skipping**, not protease cleavage:

- [YES] **2A skipping**: NPG|P junction, no protease needed
- [NO] **Protease cleavage**: Q|G or E|G, requires viral protease (excluded)

## Examples of Expected Results

### Foot-and-Mouth Disease Virus (FMDV)
- Family: Picornaviridae
- 2A motif: LLNFDLLKLAGDVESNPGP
- Downstream: P2B, P2C, P3 proteins

### Enterovirus
- Family: Picornaviridae  
- 2A motif: DIESTNPGP (highly conserved)
- Downstream: Structural proteins

## Troubleshooting

### "No sequences downloaded"
- Check your internet connection
- UniProt API may be temporarily unavailable
- Try again later

### "Low number of results"
- Increase batch_size in the script
- The default (500) is conservative for speed

### Want more results?
- Edit line: `finder.run_analysis(batch_size=1000)`
- Or higher (2000, 5000) for comprehensive search

## References

- Donnelly ML, et al. (2001) "The 'cleavage' activities of foot-and-mouth disease virus 2A site-directed mutants and naturally occurring '2A-like' sequences." J Gen Virol.
- Ryan MD, Drew J. (1994) "Foot-and-mouth disease virus 2A oligopeptide mediated cleavage of an artificial polyprotein." EMBO J.

## License

MIT License - Feel free to use and modify for research purposes.

## Author

Created for viral protein bioinformatics research.
