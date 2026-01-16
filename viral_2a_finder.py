#!/usr/bin/env python3
"""
Viral 2A Motif Finder
Identifies viral proteins with functional 2A motifs that undergo ribosomal skipping
"""

import re
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import pandas as pd
from datetime import datetime
import time

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

class Viral2AFinder:
    def __init__(self):
        self.results = []
        self.uniprot_base = "https://rest.uniprot.org/uniprotkb"
        
    def download_viral_polyproteins(self, batch_size=500):
        """
        Download viral polyproteins with chain annotations from UniProt
        """
        print("=" * 80)
        print("STEP 1: Downloading viral polyproteins from UniProt...")
        print("=" * 80)
        
        # UniProt query for viral polyproteins with chain annotations
        # Updated query syntax for new UniProt API
        query = "(taxonomy_id:10239) AND (cc_function:polyprotein)"
        
        url = f"{self.uniprot_base}/search"
        params = {
            'query': query,
            'format': 'fasta',
            'size': batch_size
        }
        
        print(f"Query: {query}")
        print(f"Fetching up to {batch_size} sequences...")
        
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            # Parse FASTA sequences
            fasta_io = StringIO(response.text)
            sequences = list(SeqIO.parse(fasta_io, "fasta"))
            
            print(f"[OK] Downloaded {len(sequences)} viral polyprotein sequences")
            return sequences
            
        except Exception as e:
            print(f"[ERROR] Error downloading sequences: {e}")
            print("Trying alternative query...")
            
            # Try alternative query
            query2 = "(taxonomy_id:10239) AND (keyword:polyprotein)"
            params['query'] = query2
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
        Returns list of (motif_sequence, position, pattern_type)
        """
        motifs = []
        
        # Pattern 1: Canonical 2A motif (most common)
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
            # Avoid duplicates from pattern1
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
        # Position after the G in NPGP (where cleavage occurs)
        cleavage_pos = motif_end_pos - 1  # -1 because we want after G, before P
        
        # Extract downstream sequence
        downstream = str(sequence)[cleavage_pos:]
        
        # Check if it starts with P (Proline)
        if downstream and downstream[0] == 'P':
            return downstream
        return None
    
    def is_known_2a_family(self, description):
        """
        Check if the virus belongs to a known 2A-containing family
        """
        description_lower = description.lower()
        
        # Check each family and its associated virus names
        for family, virus_patterns in KNOWN_2A_FAMILIES.items():
            for pattern in virus_patterns:
                if pattern.lower() in description_lower:
                    return family
        
        return None
    
    def extract_virus_name(self, description):
        """
        Extract virus name from UniProt description
        """
        # UniProt format: >sp|ACCESSION|NAME Description OS=Organism
        if 'OS=' in description:
            parts = description.split('OS=')
            if len(parts) > 1:
                organism = parts[1].split('OX=')[0].strip()
                return organism
        
        # Fallback: return first part after identifier
        parts = description.split()
        if len(parts) > 1:
            return ' '.join(parts[1:5])  # Take first few words
        return description
    
    def analyze_sequences(self, sequences):
        """
        Main analysis pipeline
        """
        print("\n" + "=" * 80)
        print("STEP 2-4: Analyzing sequences for 2A motifs...")
        print("=" * 80)
        
        total_2a_found = 0
        valid_2a_proteins = 0
        
        for idx, record in enumerate(sequences, 1):
            if idx % 50 == 0:
                print(f"  Processed {idx}/{len(sequences)} sequences...")
            
            # Find 2A motifs
            motifs = self.find_2a_motifs(record.seq)
            
            if motifs:
                total_2a_found += len(motifs)
                
                # Check if from known 2A family
                family = self.is_known_2a_family(record.description)
                virus_name = self.extract_virus_name(record.description)
                
                # Process each motif
                for motif in motifs:
                    # Extract downstream protein
                    downstream = self.extract_downstream_protein(
                        record.seq, 
                        motif['end_position']
                    )
                    
                    if downstream:
                        valid_2a_proteins += 1
                        
                        # Store result
                        result = {
                            'virus_name': virus_name,
                            'accession': record.id.split('|')[1] if '|' in record.id else record.id,
                            'family': family if family else 'Unknown',
                            '2a_motif': motif['sequence'],
                            '2a_pattern': motif['pattern'],
                            '2a_position': motif['position'],
                            'downstream_length': len(downstream),
                            'downstream_sequence': downstream,
                            'confidence': 'HIGH' if family else 'LOW',
                            'full_description': record.description
                        }
                        
                        self.results.append(result)
        
        print(f"\n[OK] Found {total_2a_found} total 2A motifs")
        print(f"[OK] Found {valid_2a_proteins} valid 2A proteins with P-terminal downstream products")
        print(f"[OK] High confidence (known families): {sum(1 for r in self.results if r['confidence'] == 'HIGH')}")
        
    def generate_outputs(self):
        """
        Generate output files with results
        """
        print("\n" + "=" * 80)
        print("STEP 5: Generating output files...")
        print("=" * 80)
        
        if not self.results:
            print("No results to save.")
            return
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 1. Summary CSV
        df = pd.DataFrame(self.results)
        csv_file = f"2a_proteins_summary_{timestamp}.csv"
        df_summary = df[['virus_name', 'accession', 'family', '2a_motif', 
                        '2a_position', 'downstream_length', 'confidence']]
        df_summary.to_csv(csv_file, index=False)
        print(f"[OK] Saved summary to: {csv_file}")
        
        # 2. Detailed Excel with all info
        excel_file = f"2a_proteins_detailed_{timestamp}.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='All Results', index=False)
            
            # High confidence only
            df_high = df[df['confidence'] == 'HIGH']
            if not df_high.empty:
                df_high.to_excel(writer, sheet_name='High Confidence', index=False)
            
            # Group by family
            df_by_family = df.groupby('family').size().reset_index(name='count')
            df_by_family.to_excel(writer, sheet_name='By Family', index=False)
        
        print(f"[OK] Saved detailed results to: {excel_file}")
        
        # 3. FASTA file with downstream proteins
        fasta_file = f"2a_downstream_proteins_{timestamp}.fasta"
        records = []
        for result in self.results:
            seq_id = f"{result['accession']}|{result['2a_motif']}|pos{result['2a_position']}"
            description = f"{result['virus_name']} | {result['family']} | {result['confidence']}"
            
            record = SeqRecord(
                Seq(result['downstream_sequence']),
                id=seq_id,
                description=description
            )
            records.append(record)
        
        SeqIO.write(records, fasta_file, "fasta")
        print(f"[OK] Saved downstream proteins to: {fasta_file}")
        
        # 4. High confidence FASTA only
        high_conf_fasta = f"2a_downstream_HIGH_CONFIDENCE_{timestamp}.fasta"
        high_conf_records = [r for i, r in enumerate(records) 
                            if self.results[i]['confidence'] == 'HIGH']
        if high_conf_records:
            SeqIO.write(high_conf_records, high_conf_fasta, "fasta")
            print(f"[OK] Saved high-confidence proteins to: {high_conf_fasta}")
        
        # 5. Summary statistics
        print("\n" + "=" * 80)
        print("SUMMARY STATISTICS")
        print("=" * 80)
        print(f"Total 2A proteins found: {len(self.results)}")
        print(f"High confidence: {sum(1 for r in self.results if r['confidence'] == 'HIGH')}")
        print(f"Low confidence: {sum(1 for r in self.results if r['confidence'] == 'LOW')}")
        print("\nBy viral family:")
        family_counts = df.groupby('family').size().sort_values(ascending=False)
        for family, count in family_counts.items():
            print(f"  {family}: {count}")
        
        print("\n2A motif patterns found:")
        motif_counts = df.groupby('2a_motif').size().sort_values(ascending=False)
        for motif, count in motif_counts.head(10).items():
            print(f"  {motif}: {count}")
    
    def run_analysis(self, batch_size=500):
        """
        Run the complete analysis pipeline
        """
        print("\n")
        print("+" + "=" * 78 + "+")
        print("|" + " " * 20 + "VIRAL 2A PROTEIN FINDER" + " " * 35 + "|")
        print("+" + "=" * 78 + "+")
        print(f"\nStarted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        # Step 1: Download sequences
        sequences = self.download_viral_polyproteins(batch_size=batch_size)
        
        if not sequences:
            print("No sequences downloaded. Exiting.")
            return
        
        # Step 2-4: Analyze sequences
        self.analyze_sequences(sequences)
        
        # Step 5: Generate outputs
        self.generate_outputs()
        
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("\nFiles generated:")
        print("  - *_summary_*.csv (overview of all results)")
        print("  - *_detailed_*.xlsx (comprehensive data with multiple sheets)")
        print("  - *_downstream_proteins_*.fasta (all downstream sequences)")
        print("  - *_HIGH_CONFIDENCE_*.fasta (curated high-confidence hits)")
        print("\n")


def main():
    """
    Main entry point
    """
    finder = Viral2AFinder()
    
    # Run analysis with configurable batch size
    # Increase batch_size for more comprehensive results
    finder.run_analysis(batch_size=500)


if __name__ == "__main__":
    main()
