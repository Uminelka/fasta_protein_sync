import time
import json
import threading
import queue
from collections import Counter
from Bio import SeqIO

class ProteinAnalyzerThreaded:
    def __init__(self, file, num_threads=4):
        self.file = file
        self.num_threads = num_threads
        
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
        
        self.task_queue = queue.Queue(maxsize=2000)
        self.data_lock = threading.Lock()

    def _worker(self):
        while True:
            item = self.task_queue.get()
            if item is None:
                self.task_queue.task_done()
                break
                
            seq_id, seq_str = item
            counts = Counter(seq_str)
            
            current_seq_results = {
                'id': seq_id,
                'hydrophobic': 0, 'hydrophilic_neutral': 0,
                'hydrophilic_positive': 0, 'hydrophilic_negative': 0,
                'not_standart': 0
            }
            
            for amk, count in counts.items():
                defined = False
                for class_amk, set_amk in self.properties.items():
                    if amk in set_amk:
                        current_seq_results[class_amk] += count
                        defined = True
                        break
                if not defined:
                    current_seq_results['not_standart'] += count

            with self.data_lock:
                self.sequences_data.append(current_seq_results)
                for key in self.results.keys():
                    self.results[key] += current_seq_results[key]
                    
            self.task_queue.task_done()

    def analyze(self):
        start_time = time.time()
        
        threads = []
        for _ in range(self.num_threads):
            t = threading.Thread(target=self._worker)
            t.start()
            threads.append(t)
            
        with open(self.file) as f:
            for sequence in SeqIO.parse(f, 'fasta'):
                seq_str = str(sequence.seq).upper()
                self.task_queue.put((sequence.id, seq_str))
                
        for _ in range(self.num_threads):
            self.task_queue.put(None)
            
        for t in threads:
            t.join()
            
        self.full_time = time.time() - start_time

    def save_result(self, out_file):
        out_data = {
            'execution_time_seconds': self.full_time,
            'total_statistics': self.results,
            'sequences': self.sequences_data
        }
        with open(out_file, 'w', encoding='utf-8') as f:
            json.dump(out_data, f, indent=4)


if __name__ == '__main__':
    test_file = 'example_10000.fasta'
    
    analyzer_th = ProteinAnalyzerThreaded(test_file, num_threads=4)
    analyzer_th.analyze()
    analyzer_th.save_result('threaded_result.json')