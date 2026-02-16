import time
import json
from collections import Counter
from Bio import SeqIO

class ProteinAnalyzer:
    def __init__(self, file):
        self.file = file
        
        self.properties = {
            'hydrophobic': set('AVILMFWCPG'),
            'hydrophilic_neutral': set('STNQY'),
            'hydrophilic_positive': set('KRH'),
            'hydrophilic_negative': set('DE'),
        }
        
        self.results = {
            'hydrophobic': 0,
            'hydrophilic_neutral': 0,
            'hydrophilic_positive': 0,
            'hydrophilic_negative': 0,
            'not_standart': 0
        }
        
        self.sequences_data = []
        self.full_time = 0
        
    def analyze(self):
        start_time = time.time()
        
        with open(self.file) as f:
            for sequence in SeqIO.parse(f, 'fasta'):
                counts = Counter(str(sequence.seq).upper())
                
                current_seq_results = {
                    'id': sequence.id,
                    'hydrophobic': 0,
                    'hydrophilic_neutral': 0,
                    'hydrophilic_positive': 0,
                    'hydrophilic_negative': 0,
                    'not_standart': 0
                }
                
                for amk, count in counts.items():
                    defined = False
                    for class_amk, set_amk in self.properties.items():
                        if amk in set_amk:
                            current_seq_results[class_amk] += count
                            self.results[class_amk] += count
                            defined = True
                            break
                    
                    if not defined:
                        current_seq_results['not_standart'] += count
                        self.results['not_standart'] += count
                
                self.sequences_data.append(current_seq_results)
                        
        end_time = time.time()
        self.full_time = end_time - start_time
        
    def save_result(self, out_file):
        out_data = {
            'execution_time_seconds': self.full_time,
            'total_statistics': self.results,
            'sequences': self.sequences_data
        }
        
        with open(out_file, 'w') as f:
            json.dump(out_data, f, indent=4)
    
if __name__ == '__main__':
    test_file = 'example_100.fasta'

    try:
        analyzer = ProteinAnalyzer(test_file)
        analyzer.analyze()
        analyzer.save_result('fasta_analyzer.json')
        
    except FileNotFoundError as e:
        print(e)