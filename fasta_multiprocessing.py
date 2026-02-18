import time
import json
import multiprocessing
from collections import Counter
from Bio import SeqIO
from functools import partial

 
def pool_worker(seq_data, properties):
    seq_id, seq_str = seq_data
    counts = Counter(seq_str)
    
    current_seq_results = {
        'id': seq_id,
        'hydrophobic': 0, 'hydrophilic_neutral': 0,
        'hydrophilic_positive': 0, 'hydrophilic_negative': 0,
        'not_standart': 0
    }
    
    for amk, count in counts.items():
        defined = False
        for class_amk, set_amk in properties.items():
            if amk in set_amk:
                current_seq_results[class_amk] += count
                defined = True
                break
        if not defined:
            current_seq_results['not_standart'] += count
    return current_seq_results

class ProteinAnalyzerPool:
    def __init__(self, file):
        self.file = file
        self.num_processes = 12
        self.properties = {
            'hydrophobic': set('AVILMFWCPG'),
            'hydrophilic_neutral': set('STNQY'),
            'hydrophilic_positive': set('KRH'),
            'hydrophilic_negative': set('DE'),
        }
        self.results = {key: 0 for key in list(self.properties.keys()) + ['not_standart']}
        self.sequences_data = []
        self.full_time = 0

    def analyze(self):
        start_time = time.time()

         
        def sequence_generator():
            with open(self.file) as f:
                for sequence in SeqIO.parse(f, 'fasta'):
                    yield (sequence.id, str(sequence.seq).upper())

         
        with multiprocessing.Pool(processes=self.num_processes) as pool:
             
            worker_func = partial(pool_worker, properties=self.properties)
            
            
            for res in pool.imap_unordered(worker_func, sequence_generator(), chunksize=50):
                self.sequences_data.append(res)
                for key in self.results:
                    self.results[key] += res[key]

        self.full_time = time.time() - start_time

    def save_result(self, out_file):
        out_data = {
            'execution_time_seconds': self.full_time,
            'total_statistics': self.results,
            'sequences': self.sequences_data
        }
        with open(out_file, 'w', encoding='utf-8') as f:
            json.dump(out_data, f, indent=4)

if __name__ == "__main__":
    test_file = "example_100000.fasta"
    analyzer = ProteinAnalyzerPool(test_file)
    analyzer.analyze()
    analyzer.save_result("pool_result.json")
