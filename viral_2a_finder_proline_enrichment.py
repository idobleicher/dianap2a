#!/usr/bin/env python3
"""
Viral 2A Motif Finder - Proline Position Analysis
Filters for P at position 1 and D/E/T at position 2
Includes enrichment analysis of residues following Proline
"""

import re
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import pandas as pd
import numpy as np
from datetime import datetime
import time
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict

# Known 2A-containing viral families (high confidence)
KNOWN_2A_FAMILIES = {
    'Picornaviridae': [
        'foot-and-mouth disease virus',
        'fmdv',
        'enterovirus',
        'poliovirus',
        'rhinovirus',
        'hepatovirus',
        'cardiovirus',
        'aphthovirus',
        'kobuvirus',
        'teschovirus',
        'seneca valley virus',
        'ljunganvirus'
    ],
    'Dicistroviridae': [
        'cricket paralysis virus',
        'triatoma virus',
        'drosophila c virus'
    ],
    'Iflaviridae': [
        'iflavirus',
        'deformed wing virus',
        'sacbrood virus'
    ],
    'Unclassified_picorna-like': [
        'bee paralysis virus',
        'picorna-like virus'
    ]
}

# Standard amino acids
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

class Viral2AFinder:
    def __init__(self, analysis_window=20):
        """
        Initialize finder with configurable analysis window
        analysis_window: number of positions to analyze after P
        """
        self.results = []
        self.all_downstream = []  # Store all downstream sequences for enrichment
        self.all_polyproteins = []  # Store ALL polyproteins for virome-wide analysis
        self.uniprot_base = "https://rest.uniprot.org/uniprotkb"
        self.analysis_window = analysis_window
        
    def download_viral_polyproteins(self, total_wanted=2000):
        """
        Download viral polyproteins with chain annotations from UniProt
        Uses pagination to get more than 500 sequences
        """
        print("=" * 80)
        print("STEP 1: Downloading viral polyproteins from UniProt...")
        print("=" * 80)
        
        query = "(taxonomy_id:10239) AND (cc_function:polyprotein)"
        url = f"{self.uniprot_base}/search"
        
        all_sequences = []
        size_per_request = 500  # Max safe size per request
        cursor = None
        
        print(f"Query: {query}")
        print(f"Fetching sequences in batches of {size_per_request}...")
        
        try:
            while len(all_sequences) < total_wanted:
                params = {
                    'query': query,
                    'format': 'fasta',
                    'size': size_per_request
                }
                
                if cursor:
                    params['cursor'] = cursor
                
                response = requests.get(url, params=params)
                response.raise_for_status()
                
                # Parse FASTA sequences
                fasta_io = StringIO(response.text)
                batch_sequences = list(SeqIO.parse(fasta_io, "fasta"))
                
                if not batch_sequences:
                    break  # No more sequences
                
                all_sequences.extend(batch_sequences)
                print(f"  Downloaded {len(all_sequences)} sequences so far...")
                
                # Check if there's a cursor for next page
                if 'Link' in response.headers:
                    link_header = response.headers['Link']
                    if 'cursor=' in link_header:
                        # Extract cursor from Link header
                        cursor_match = re.search(r'cursor=([^&>]+)', link_header)
                        if cursor_match:
                            cursor = cursor_match.group(1)
                        else:
                            break
                    else:
                        break
                else:
                    break
                
                # Small delay to be nice to the API
                time.sleep(0.5)
            
            print(f"[OK] Downloaded {len(all_sequences)} viral polyprotein sequences")
            return all_sequences
            
        except Exception as e:
            print(f"[ERROR] Error downloading sequences: {e}")
            if all_sequences:
                print(f"[OK] Returning {len(all_sequences)} sequences downloaded so far")
                return all_sequences
            
            print("Trying alternative query...")
            
            # Try alternative query with single batch
            query2 = "(taxonomy_id:10239) AND (keyword:polyprotein)"
            params = {'query': query2, 'format': 'fasta', 'size': 500}
            try:
                response = requests.get(url, params=params)
                response.raise_for_status()
                fasta_io = StringIO(response.text)
                sequences = list(SeqIO.parse(fasta_io, "fasta"))
                print(f"[OK] Downloaded {len(sequences)} sequences with alternative query")
                return sequences
            except:
                print("[ERROR] Alternative query also failed")
                return []
    
    def find_2a_motifs(self, sequence):
        """
        Find 2A motifs in a protein sequence
        """
        motifs = []
        
        # Pattern 1: Canonical 2A motif
        pattern1 = re.compile(r'D[VIL]E[ST]NPGP')
        for match in pattern1.finditer(str(sequence)):
            motifs.append({
                'sequence': match.group(),
                'position': match.start(),
                'end_position': match.end(),
                'pattern': 'D[VIL]E[ST]NPGP (canonical)'
            })
        
        # Pattern 2: Broader 2A motif
        pattern2 = re.compile(r'D[A-Z][A-Z]E[A-Z]NPGP')
        for match in pattern2.finditer(str(sequence)):
            if not any(m['position'] == match.start() for m in motifs):
                motifs.append({
                    'sequence': match.group(),
                    'position': match.start(),
                    'end_position': match.end(),
                    'pattern': 'D..E.NPGP (broader)'
                })
        
        return motifs
    
    def extract_downstream_protein(self, sequence, motif_end_pos):
        """
        Extract downstream protein sequence after NPG|P junction
        The cleavage occurs between G and P
        """
        cleavage_pos = motif_end_pos - 1
        downstream = str(sequence)[cleavage_pos:]
        
        if downstream and downstream[0] == 'P':
            return downstream
        return None
    
    def check_position_2(self, downstream_seq):
        """
        Check if position 2 (immediately after P) is D/E/T
        
        Returns: match info if position 2 is D/E/T, None otherwise
        """
        if len(downstream_seq) < 2:
            return None
            
        # Position 1 should be P (already validated)
        # Position 2 is index 1
        if downstream_seq[1] in ['D', 'E', 'T']:
            # Get context (first 5-10 residues)
            context_length = min(10, len(downstream_seq))
            context = downstream_seq[:context_length]
            
            return {
                'position_1': downstream_seq[0],  # Should be P
                'position_2': downstream_seq[1],  # D/E/T
                'dipeptide': downstream_seq[:2],
                'context': context
            }
        
        return None
    
    def is_known_2a_family(self, description):
        """
        Check if the virus belongs to a known 2A-containing family
        """
        description_lower = description.lower()
        
        for family, virus_patterns in KNOWN_2A_FAMILIES.items():
            for pattern in virus_patterns:
                if pattern.lower() in description_lower:
                    return family
        
        return None
    
    def extract_virus_name(self, description):
        """
        Extract virus name from UniProt description
        """
        if 'OS=' in description:
            parts = description.split('OS=')
            if len(parts) > 1:
                organism = parts[1].split('OX=')[0].strip()
                return organism
        
        parts = description.split()
        if len(parts) > 1:
            return ' '.join(parts[1:5])
        return description
    
    def analyze_sequences(self, sequences):
        """
        Main analysis pipeline
        """
        print("\n" + "=" * 80)
        print("STEP 2-4: Analyzing sequences for 2A motifs...")
        print("=" * 80)
        print("Filter: P at position 1, D/E/T at position 2")
        print(f"Storing ALL {len(sequences)} polyproteins for virome-wide enrichment analysis")
        
        total_2a_found = 0
        valid_2a_proteins = 0
        filtered_proteins = 0
        
        for idx, record in enumerate(sequences, 1):
            if idx % 50 == 0:
                print(f"  Processed {idx}/{len(sequences)} sequences...")
            
            # Store ALL polyproteins for enrichment analysis
            family = self.is_known_2a_family(record.description)
            virus_name = self.extract_virus_name(record.description)
            self.all_polyproteins.append({
                'sequence': str(record.seq),
                'family': family if family else 'Unknown',
                'virus_name': virus_name
            })
            
            motifs = self.find_2a_motifs(record.seq)
            
            if motifs:
                total_2a_found += len(motifs)
                
                family = self.is_known_2a_family(record.description)
                virus_name = self.extract_virus_name(record.description)
                
                for motif in motifs:
                    downstream = self.extract_downstream_protein(
                        record.seq, 
                        motif['end_position']
                    )
                    
                    if downstream:
                        valid_2a_proteins += 1
                        
                        # Store ALL downstream sequences for enrichment analysis
                        self.all_downstream.append({
                            'sequence': downstream,
                            'family': family if family else 'Unknown',
                            'virus_name': virus_name,
                            'motif_start': motif['position'],
                            'protein_sequence': str(record.seq)
                        })
                        
                        # Check position 2 filter
                        position_match = self.check_position_2(downstream)
                        
                        if position_match:
                            filtered_proteins += 1
                            
                            result = {
                                'virus_name': virus_name,
                                'accession': record.id.split('|')[1] if '|' in record.id else record.id,
                                'family': family if family else 'Unknown',
                                '2a_motif': motif['sequence'],
                                '2a_pattern': motif['pattern'],
                                '2a_position': motif['position'],
                                'downstream_length': len(downstream),
                                'downstream_sequence': downstream,
                                # Filtered fields
                                'position_1': position_match['position_1'],
                                'position_2': position_match['position_2'],
                                'dipeptide': position_match['dipeptide'],
                                'context': position_match['context'],
                                'confidence': 'HIGH' if family else 'LOW',
                                'full_description': record.description
                            }
                            
                            self.results.append(result)
        
        print(f"\n[OK] Found {total_2a_found} total 2A motifs")
        print(f"[OK] Found {valid_2a_proteins} valid 2A proteins with P-terminal downstream products")
        print(f"[OK] FILTERED: {filtered_proteins} proteins with P at pos1 and D/E/T at pos2")
        print(f"[OK] High confidence (known families): {sum(1 for r in self.results if r['confidence'] == 'HIGH')}")
        
    def calculate_enrichment(self):
        """
        Calculate amino acid enrichment at each position in 2A downstream proteins
        Compare to background amino acid frequency
        """
        print("\n" + "=" * 80)
        print("ENRICHMENT ANALYSIS: Amino acid frequencies in 2A downstream proteins")
        print("=" * 80)
        print(f"Analyzing {len(self.all_downstream)} 2A downstream protein sequences")
        
        if not self.all_downstream:
            print("No sequences to analyze")
            return None, None, None, None
        
        # Calculate position-specific frequencies in 2A downstream proteins
        position_counts = defaultdict(lambda: Counter())
        
        for entry in self.all_downstream:
            seq = entry['sequence']
            # Analyze first N positions
            for i in range(min(len(seq), self.analysis_window)):
                aa = seq[i]
                if aa in AMINO_ACIDS:  # Skip non-standard amino acids
                    position_counts[i+1][aa] += 1  # 1-indexed positions
        
        # Calculate background frequency (all positions combined in 2A downstream proteins)
        background = Counter()
        for entry in self.all_downstream:
            for aa in entry['sequence']:  # Use entire sequence for background
                if aa in AMINO_ACIDS:
                    background[aa] += 1
        
        total_background = sum(background.values())
        background_freq = {aa: background[aa] / total_background for aa in AMINO_ACIDS}
        
        # Create frequency matrix
        positions = sorted(position_counts.keys())
        freq_matrix = np.zeros((len(AMINO_ACIDS), len(positions)))
        enrichment_matrix = np.zeros((len(AMINO_ACIDS), len(positions)))
        
        for pos_idx, pos in enumerate(positions):
            total = sum(position_counts[pos].values())
            for aa_idx, aa in enumerate(AMINO_ACIDS):
                count = position_counts[pos][aa]
                freq = count / total if total > 0 else 0
                freq_matrix[aa_idx, pos_idx] = freq
                
                # Calculate enrichment (log2 fold change vs background)
                if background_freq[aa] > 0:
                    enrichment = freq / background_freq[aa]
                    enrichment_matrix[aa_idx, pos_idx] = np.log2(enrichment + 0.01)  # Add pseudocount
                else:
                    enrichment_matrix[aa_idx, pos_idx] = 0
        
        return freq_matrix, enrichment_matrix, positions, background_freq
    
    def plot_enrichment_heatmap(self, enrichment_matrix, positions, output_file):
        """
        Create heatmap of amino acid enrichment at each position
        """
        print("\nGenerating enrichment heatmap...")
        
        plt.figure(figsize=(14, 10))
        
        # Create heatmap with soft warm colors (red, orange, pink, beige)
        import matplotlib.colors as mcolors
        
        # Custom colormap: beige (low) -> peach/pink (medium) -> coral/red (high)
        colors = ['#F5DEB3', '#FFE4B5', '#FFDAB9', '#FFC0CB', '#FFB6C1', 
                  '#FF9999', '#FF7F7F', '#FF6B6B', '#E74C3C', '#C0392B']
        n_bins = 100
        cmap = mcolors.LinearSegmentedColormap.from_list('warm_soft', colors, N=n_bins)
        
        # Custom colormap: very light purple/lavender (depleted) -> beige (neutral) -> coral/red (enriched)
        colors = ['#F8F0FF', '#F0E6FF', '#E6D9F0', '#E8E0F0', '#F5F5DC', '#FFEFD5', 
                  '#FFE4B5', '#FFDAB9', '#FFC0CB', '#FFB6C1', '#FF9999', 
                  '#FF7F7F', '#FF6B6B', '#E74C3C', '#C0392B']
        cmap = mcolors.LinearSegmentedColormap.from_list('lavender_beige_red', colors, N=100)
        
        ax = sns.heatmap(
            enrichment_matrix,
            xticklabels=positions,
            yticklabels=AMINO_ACIDS,
            cmap=cmap,  # Custom: light blue -> beige -> red
            center=0,
            vmin=-3,
            vmax=3,
            cbar_kws={'label': 'Log2 Enrichment vs Background'},
            linewidths=0.2,
            linecolor='#F5F5DC',  # Beige lines
            annot=False
        )
        
        plt.title(f'Amino Acid Enrichment in 2A Downstream Proteins (n={len(self.all_downstream)})\n(First {self.analysis_window} positions after 2A cleavage)', 
                  fontsize=14, fontweight='bold', pad=20, color='#8B4513')
        plt.xlabel('Position', fontsize=12, fontweight='bold', color='#A0522D')
        plt.ylabel('Amino Acid', fontsize=12, fontweight='bold', color='#A0522D')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved enrichment heatmap to: {output_file}")
        
    def plot_frequency_heatmap(self, freq_matrix, positions, output_file):
        """
        Create heatmap of raw amino acid frequencies
        """
        print("Generating frequency heatmap...")
        
        plt.figure(figsize=(14, 10))
        
        # Custom warm colormap: beige -> peach -> coral -> pink -> red/orange
        import matplotlib.colors as mcolors
        colors = ['#FFF8DC', '#FFEBCD', '#FFE4C4', '#FFDAB9', '#FFDEAD',
                  '#FFB6C1', '#FFA07A', '#FF8C69', '#FF7F50', '#FF6347', 
                  '#FF4500', '#DC143C']
        cmap = mcolors.LinearSegmentedColormap.from_list('warm_gradient', colors, N=100)
        
        # Start from position 2 (skip position 1 which is always P)
        freq_matrix_from_2 = freq_matrix[:, 1:]  # Skip first column
        positions_from_2 = positions[1:]  # Skip first position
        
        ax = sns.heatmap(
            freq_matrix_from_2 * 100,  # Convert to percentages
            xticklabels=positions_from_2,
            yticklabels=AMINO_ACIDS,
            cmap=cmap,  # Custom warm gradient
            cbar_kws={'label': 'Frequency (%)'},
            linewidths=0.2,
            linecolor='#FFF8DC',  # Cornsilk/beige lines
            annot=False
        )
        
        plt.title(f'Amino Acid Frequency in 2A Downstream Proteins (n={len(self.all_downstream)})\n(Positions 2-{self.analysis_window} after 2A cleavage)', 
                  fontsize=14, fontweight='bold', pad=20, color='#8B4513')
        plt.xlabel('Position', fontsize=12, fontweight='bold', color='#A0522D')
        plt.ylabel('Amino Acid', fontsize=12, fontweight='bold', color='#A0522D')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved frequency heatmap to: {output_file}")
    
    def generate_enrichment_table(self, freq_matrix, enrichment_matrix, positions, background_freq, timestamp):
        """
        Generate detailed enrichment table
        """
        # Create DataFrame for each position
        data = []
        for pos_idx, pos in enumerate(positions):
            for aa_idx, aa in enumerate(AMINO_ACIDS):
                data.append({
                    'Position': pos,
                    'Amino_Acid': aa,
                    'Frequency': freq_matrix[aa_idx, pos_idx],
                    'Frequency_Percent': freq_matrix[aa_idx, pos_idx] * 100,
                    'Background_Freq': background_freq[aa],
                    'Log2_Enrichment': enrichment_matrix[aa_idx, pos_idx],
                    'Fold_Change': 2 ** enrichment_matrix[aa_idx, pos_idx]
                })
        
        df_enrichment = pd.DataFrame(data)
        
        # Save to Excel
        excel_file = f"2a_enrichment_analysis_{timestamp}.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df_enrichment.to_excel(writer, sheet_name='All Positions', index=False)
            
            # Position 2 only (our filter position)
            df_pos2 = df_enrichment[df_enrichment['Position'] == 2].sort_values('Frequency', ascending=False)
            df_pos2.to_excel(writer, sheet_name='Position 2 (Filter)', index=False)
            
            # Top enriched at each position
            df_top = df_enrichment.sort_values('Log2_Enrichment', ascending=False).head(50)
            df_top.to_excel(writer, sheet_name='Top Enriched', index=False)
        
        print(f"[OK] Saved enrichment analysis to: {excel_file}")
        
        return df_enrichment
    
    def analyze_upstream_phenylalanine(self, timestamp):
        """
        Analyze F (Phenylalanine) in upstream protein (N-terminal region)
        """
        print("\n" + "=" * 80)
        print("UPSTREAM PROTEIN F ANALYSIS (N-terminal region)")
        print("=" * 80)
        
        if not self.all_downstream:
            print("No sequences to analyze")
            return None, None
        
        # Extract upstream sequences (N-terminal region, first 20 positions)
        window = 20
        f_counts_upstream = defaultdict(int)
        total_counts_upstream = defaultdict(int)
        
        for entry in self.all_downstream:
            motif_start = entry['motif_start']
            full_seq = entry['protein_sequence']
            
            # Get upstream sequence (entire sequence before 2A motif)
            upstream_seq = full_seq[:motif_start]
            
            # Analyze N-terminal region (first 20 positions of upstream protein)
            for i in range(min(len(upstream_seq), window)):
                pos = i + 1  # 1, 2, 3, etc. (1-indexed)
                aa = upstream_seq[i]
                total_counts_upstream[pos] += 1
                if aa == 'F':
                    f_counts_upstream[pos] += 1
        
        # Calculate frequencies
        f_freq_upstream = {}
        for pos in range(1, window + 1):
            if total_counts_upstream[pos] > 0:
                f_freq_upstream[pos] = (f_counts_upstream[pos] / total_counts_upstream[pos]) * 100
            else:
                f_freq_upstream[pos] = 0
        
        # Calculate background
        total_f = sum(f_counts_upstream.values())
        total_aa = sum(total_counts_upstream.values())
        background_f_upstream = (total_f / total_aa) * 100 if total_aa > 0 else 0
        
        print(f"\nUpstream Background F frequency (N-terminal region): {background_f_upstream:.2f}%")
        print(f"\nF frequency at each upstream position (N-terminal, positions 1-{window}):")
        print("-" * 70)
        print(f"{'Position':<12} {'F Count':<12} {'Total':<12} {'F %':<12} {'vs Background':<15}")
        print("-" * 70)
        
        for pos in range(1, window + 1):
            if pos in f_freq_upstream:
                f_pct = f_freq_upstream[pos]
                fold_change = f_pct / background_f_upstream if background_f_upstream > 0 else 0
                enrichment = "+" if fold_change > 1.5 else ("" if fold_change > 0.5 else "-")
                print(f"{pos:<12} {f_counts_upstream[pos]:<12} {total_counts_upstream[pos]:<12} {f_pct:<12.2f} {fold_change:.2f}x {enrichment}")
        
        return f_freq_upstream, f_counts_upstream, total_counts_upstream, background_f_upstream
    
    def analyze_phenylalanine_distribution(self, timestamp):
        """
        Analyze F (Phenylalanine) distribution across all positions
        """
        print("\n" + "=" * 80)
        print("PHENYLALANINE (F) DISTRIBUTION ANALYSIS")
        print("=" * 80)
        
        if not self.all_downstream:
            print("No sequences to analyze")
            return
        
        # Count F at each position
        f_counts = defaultdict(int)
        total_counts = defaultdict(int)
        
        for entry in self.all_downstream:
            seq = entry['sequence']
            for i in range(min(len(seq), self.analysis_window)):
                aa = seq[i]
                total_counts[i+1] += 1
                if aa == 'F':
                    f_counts[i+1] += 1
        
        # Calculate frequencies
        f_freq_by_position = {}
        for pos in range(1, self.analysis_window + 1):
            if total_counts[pos] > 0:
                f_freq_by_position[pos] = (f_counts[pos] / total_counts[pos]) * 100
            else:
                f_freq_by_position[pos] = 0
        
        # Calculate background F frequency (all positions)
        total_f = sum(f_counts.values())
        total_aa = sum(total_counts.values())
        background_f_freq = (total_f / total_aa) * 100 if total_aa > 0 else 0
        
        # Print results
        print(f"\nBackground F frequency (all positions): {background_f_freq:.2f}%")
        print(f"\nF frequency at each position:")
        print("-" * 60)
        print(f"{'Position':<10} {'F Count':<12} {'Total':<12} {'F %':<12} {'vs Background':<15}")
        print("-" * 60)
        
        for pos in range(1, self.analysis_window + 1):
            if pos in f_freq_by_position:
                f_pct = f_freq_by_position[pos]
                fold_change = f_pct / background_f_freq if background_f_freq > 0 else 0
                enrichment = "+" if fold_change > 1.5 else ("" if fold_change > 0.5 else "-")
                print(f"{pos:<10} {f_counts[pos]:<12} {total_counts[pos]:<12} {f_pct:<12.2f} {fold_change:.2f}x {enrichment}")
        
        # Create visualization
        import matplotlib.pyplot as plt
        
        positions = list(range(2, self.analysis_window + 1))  # Skip position 1
        f_frequencies = [f_freq_by_position.get(pos, 0) for pos in positions]
        
        plt.figure(figsize=(14, 6))
        bars = plt.bar(positions, f_frequencies, color='#FF6B6B', alpha=0.8, edgecolor='#C0392B', linewidth=1.5)
        plt.axhline(y=background_f_freq, color='#3498DB', linestyle='--', linewidth=2, 
                    label=f'Background F frequency ({background_f_freq:.1f}%)')
        
        # Highlight positions 2-4
        for i, pos in enumerate(positions):
            if pos in [2, 3, 4]:
                bars[i].set_color('#DC143C')
                bars[i].set_alpha(1.0)
        
        plt.xlabel('Position', fontsize=12, fontweight='bold', color='#8B4513')
        plt.ylabel('F (Phenylalanine) Frequency (%)', fontsize=12, fontweight='bold', color='#8B4513')
        plt.title(f'Phenylalanine (F) Distribution in 2A Downstream Proteins (n={len(self.all_downstream)})\nPositions 2-{self.analysis_window}',
                  fontsize=14, fontweight='bold', color='#8B4513', pad=20)
        plt.legend(fontsize=11)
        plt.grid(axis='y', alpha=0.3, linestyle=':', color='gray')
        plt.xticks(positions, [str(p) for p in positions])
        
        output_file = f"2a_phenylalanine_analysis_{timestamp}.png"
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n[OK] Saved F distribution plot to: {output_file}")
        
        # Save detailed data to CSV
        csv_data = []
        for pos in range(1, self.analysis_window + 1):
            f_pct = f_freq_by_position.get(pos, 0)
            fold_change = f_pct / background_f_freq if background_f_freq > 0 else 0
            csv_data.append({
                'Position': pos,
                'F_Count': f_counts[pos],
                'Total_Count': total_counts[pos],
                'F_Frequency_%': f_pct,
                'Background_%': background_f_freq,
                'Fold_Change': fold_change,
                'Enriched': 'Yes' if fold_change > 1.5 else 'No'
            })
        
        df_f = pd.DataFrame(csv_data)
        csv_file = f"2a_phenylalanine_analysis_{timestamp}.csv"
        df_f.to_csv(csv_file, index=False)
        print(f"[OK] Saved F analysis data to: {csv_file}")
        
        # Now analyze upstream (before 2A)
        f_freq_upstream, f_counts_upstream, total_counts_upstream, bg_upstream = self.analyze_upstream_phenylalanine(timestamp)
        
        # Create comparison plot
        print("\n" + "=" * 80)
        print("COMPARISON: Upstream vs Downstream F Distribution")
        print("=" * 80)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
        
        # Upstream plot (positions 1-20 from N-terminus)
        positions_up = list(range(1, 21))
        f_freq_up = [f_freq_upstream.get(pos, 0) for pos in positions_up]
        
        bars_up = ax1.bar(positions_up, f_freq_up, color='#4A90E2', alpha=0.8, edgecolor='#2E5C8A', linewidth=1.5)
        # Highlight position 2 in upstream
        if len(bars_up) > 1:
            bars_up[1].set_color('#2E5C8A')  # Position 2
            bars_up[1].set_alpha(1.0)
        
        ax1.axhline(y=bg_upstream, color='#E74C3C', linestyle='--', linewidth=2, 
                    label=f'Background ({bg_upstream:.1f}%)')
        ax1.set_xlabel('Position (N-terminal)', fontsize=12, fontweight='bold', color='#8B4513')
        ax1.set_ylabel('F Frequency (%)', fontsize=12, fontweight='bold', color='#8B4513')
        ax1.set_title('UPSTREAM Protein\n(N-terminal region)', fontsize=13, fontweight='bold', color='#2E5C8A')
        ax1.legend(fontsize=10)
        ax1.grid(axis='y', alpha=0.3, linestyle=':', color='gray')
        ax1.set_ylim(0, 105)
        ax1.set_xticks(positions_up)
        ax1.set_xticklabels([str(p) for p in positions_up])
        
        # Downstream plot (positions 2-20)
        positions_down = list(range(2, 21))
        f_freq_down = [f_freq_by_position.get(pos, 0) for pos in positions_down]
        
        bars = ax2.bar(positions_down, f_freq_down, color='#FF6B6B', alpha=0.8, edgecolor='#C0392B', linewidth=1.5)
        # Highlight positions 2-4 and 11
        for i, pos in enumerate(positions_down):
            if pos in [2, 3, 4, 11]:
                bars[i].set_color('#DC143C')
                bars[i].set_alpha(1.0)
        
        ax2.axhline(y=background_f_freq, color='#3498DB', linestyle='--', linewidth=2, 
                    label=f'Background ({background_f_freq:.1f}%)')
        ax2.set_xlabel('Position (after 2A)', fontsize=12, fontweight='bold', color='#8B4513')
        ax2.set_ylabel('F Frequency (%)', fontsize=12, fontweight='bold', color='#8B4513')
        ax2.set_title('DOWNSTREAM Protein\n(After 2A Cleavage)', fontsize=13, fontweight='bold', color='#C0392B')
        ax2.legend(fontsize=10)
        ax2.grid(axis='y', alpha=0.3, linestyle=':', color='gray')
        ax2.set_ylim(0, 105)
        ax2.set_xticks(positions_down)
        ax2.set_xticklabels([str(p) for p in positions_down])
        
        plt.suptitle(f'Phenylalanine (F) Distribution: Upstream vs Downstream of 2A (n={len(self.all_downstream)})',
                     fontsize=15, fontweight='bold', color='#8B4513', y=1.02)
        
        comparison_file = f"2a_phenylalanine_upstream_vs_downstream_{timestamp}.png"
        plt.tight_layout()
        plt.savefig(comparison_file, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved comparison plot to: {comparison_file}")
        
        # Print key comparison
        print("\n" + "=" * 80)
        print("KEY COMPARISON: Position 2 (upstream N-terminal) vs Position 2 (downstream N-terminal)")
        print("=" * 80)
        f_at_up2 = f_freq_upstream.get(2, 0)
        f_at_down2 = f_freq_by_position.get(2, 0)
        print(f"F at position 2 (upstream, 2nd from N-terminus): {f_at_up2:.2f}%")
        print(f"F at position 2 (downstream, 2nd from N-terminus after 2A): {f_at_down2:.2f}%")
        print(f"Difference: {f_at_down2 - f_at_up2:.2f}% {'(downstream enriched)' if f_at_down2 > f_at_up2 else '(upstream enriched)'}")
        if f_at_up2 > 0:
            print(f"Fold difference: {f_at_down2 / f_at_up2:.2f}x")
        else:
            print("Fold difference: Infinite (no F at position 2 upstream)")
        
        # Save comparison data
        comparison_data = {
            'Position_Type': ['Upstream_pos2', 'Downstream_pos2', 'Difference', 'Fold_Change'],
            'F_Frequency_%': [f_at_up2, f_at_down2, f_at_down2 - f_at_up2, 
                              f_at_down2 / f_at_up2 if f_at_up2 > 0 else float('inf')]
        }
        df_comp = pd.DataFrame(comparison_data)
        comp_csv = f"2a_phenylalanine_comparison_{timestamp}.csv"
        df_comp.to_csv(comp_csv, index=False)
        print(f"[OK] Saved comparison data to: {comp_csv}")
        
        return f_freq_by_position, background_f_freq
    
    def generate_outputs(self):
        """
        Generate output files with results
        """
        print("\n" + "=" * 80)
        print("STEP 5: Generating output files...")
        print("=" * 80)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # NEW: Phenylalanine distribution analysis
        self.analyze_phenylalanine_distribution(timestamp)
        
        # Enrichment analysis on ALL downstream proteins
        freq_matrix, enrichment_matrix, positions, background_freq = self.calculate_enrichment()
        
        if freq_matrix is not None:
            # Generate heatmaps
            self.plot_enrichment_heatmap(enrichment_matrix, positions, 
                                        f"2a_enrichment_heatmap_{timestamp}.png")
            self.plot_frequency_heatmap(freq_matrix, positions, 
                                       f"2a_frequency_heatmap_{timestamp}.png")
            
            # Generate enrichment table
            df_enrichment = self.generate_enrichment_table(freq_matrix, enrichment_matrix, 
                                                          positions, background_freq, timestamp)
        
        # Now generate outputs for FILTERED results
        if not self.results:
            print("\nNo filtered results to save.")
            print(f"Total downstream proteins analyzed for enrichment: {len(self.all_downstream)}")
            return
        
        # 1a. ALL 2A downstream proteins table
        if self.all_downstream:
            df_all_2a = pd.DataFrame([{
                'virus_name': entry['virus_name'],
                'family': entry['family'],
                'downstream_sequence': entry['sequence'],
                'downstream_length': len(entry['sequence']),
                'position_1': entry['sequence'][0] if len(entry['sequence']) > 0 else '',
                'position_2': entry['sequence'][1] if len(entry['sequence']) > 1 else '',
                'position_3': entry['sequence'][2] if len(entry['sequence']) > 2 else '',
                'first_10_residues': entry['sequence'][:10],
                'first_20_residues': entry['sequence'][:20]
            } for entry in self.all_downstream])
            
            all_2a_csv = f"2a_ALL_downstream_{timestamp}.csv"
            df_all_2a.to_csv(all_2a_csv, index=False)
            print(f"\n[OK] Saved ALL 2A downstream proteins to: {all_2a_csv}")
        
        # 1b. Summary CSV (filtered)
        df = pd.DataFrame(self.results)
        csv_file = f"2a_PD_PE_PT_filtered_{timestamp}.csv"
        df_summary = df[['virus_name', 'accession', 'family', '2a_motif', 
                        '2a_position', 'dipeptide', 'position_2', 'context',
                        'downstream_length', 'confidence']]
        df_summary.to_csv(csv_file, index=False)
        print(f"[OK] Saved P[D/E/T] filtered summary to: {csv_file}")
        
        # 2. Detailed Excel
        excel_file = f"2a_PD_PE_PT_detailed_{timestamp}.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Filtered Results', index=False)
            
            if not df.empty:
                # By dipeptide
                df_dipep = df.groupby('dipeptide').size().reset_index(name='count')
                df_dipep.to_excel(writer, sheet_name='By Dipeptide', index=False)
                
                # By family
                df_family = df.groupby('family').size().reset_index(name='count')
                df_family.to_excel(writer, sheet_name='By Family', index=False)
        
        print(f"[OK] Saved filtered detailed results to: {excel_file}")
        
        # 3. FASTA file
        fasta_file = f"2a_PD_PE_PT_proteins_{timestamp}.fasta"
        records = []
        for result in self.results:
            seq_id = f"{result['accession']}|{result['2a_motif']}|{result['dipeptide']}"
            description = f"{result['virus_name']} | {result['family']} | {result['dipeptide']} | {result['confidence']}"
            
            record = SeqRecord(
                Seq(result['downstream_sequence']),
                id=seq_id,
                description=description
            )
            records.append(record)
        
        SeqIO.write(records, fasta_file, "fasta")
        print(f"[OK] Saved filtered proteins to: {fasta_file}")
        
        # Summary statistics
        print("\n" + "=" * 80)
        print("SUMMARY STATISTICS")
        print("=" * 80)
        print(f"Total viral polyproteins analyzed: {len(self.all_polyproteins)}")
        print(f"2A downstream proteins found: {len(self.all_downstream)}")
        print(f"Proteins matching P[D/E/T] pattern: {len(self.results)}")
        print(f"High confidence: {sum(1 for r in self.results if r['confidence'] == 'HIGH')}")
        
        if not df.empty:
            print("\nDipeptide distribution (P + position 2):")
            dipep_counts = df.groupby('dipeptide').size().sort_values(ascending=False)
            for dipep, count in dipep_counts.items():
                print(f"  {dipep}: {count}")
            
            print("\nBy viral family:")
            family_counts = df.groupby('family').size().sort_values(ascending=False)
            for family, count in family_counts.items():
                print(f"  {family}: {count}")
    
    def run_analysis(self, total_sequences=2000):
        """
        Run the complete analysis pipeline
        """
        print("\n")
        print("+" + "=" * 78 + "+")
        print("|" + " " * 12 + "VIRAL 2A PROTEIN FINDER - P[D/E/T] FILTER" + " " * 23 + "|")
        print("|" + " " * 18 + "WITH ENRICHMENT ANALYSIS" + " " * 36 + "|")
        print("+" + "=" * 78 + "+")
        print(f"\nStarted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Filter: P at position 1, D/E/T at position 2")
        print(f"Enrichment analysis window: {self.analysis_window} positions\n")
        
        # Download sequences with pagination
        sequences = self.download_viral_polyproteins(total_wanted=total_sequences)
        
        if not sequences:
            print("No sequences downloaded. Exiting.")
            return
        
        # Analyze sequences
        self.analyze_sequences(sequences)
        
        # Generate outputs
        self.generate_outputs()
        
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("\nFiles generated:")
        print("  - Enrichment heatmaps (PNG)")
        print("  - Enrichment analysis table (XLSX)")
        print("  - Filtered results (CSV, XLSX, FASTA)")
        print("\n")


def main():
    """
    Main entry point
    """
    # Initialize with analysis window (positions to analyze after P)
    finder = Viral2AFinder(analysis_window=20)
    
    # Run analysis to find ALL 2A sequences (downloads up to 10000 sequences)
    finder.run_analysis(total_sequences=10000)


if __name__ == "__main__":
    main()
