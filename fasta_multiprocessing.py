import time
import json
import threading
import multiprocessing
from collections import Counter
from Bio import SeqIO

def mp_worker(task_queue, result_queue, properties):
    while True:
        item = task_queue.get()
        if item is None:
            result_queue.put(None)
            break
            
        seq_id, seq_str = item
        counts = Counter(seq_str)
        
        current_seq_results = {
            'id': seq_id,
            'hydrophobic': 0,
            'hydrophilic_neutral': 0,
            'hydrophilic_positive': 0,
            'hydrophilic_negative': 0,
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
                
        result_queue.put(current_seq_results)


class ProteinAnalyzerMultiprocess:
    def __init__(self, file):
        self.file = file
        self.num_processes = multiprocessing.cpu_count()
        
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
        
        task_queue = multiprocessing.Queue(maxsize=2000)
        result_queue = multiprocessing.Queue()
        
        processes = []
        for _ in range(self.num_processes):
            p = multiprocessing.Process(
                target=mp_worker, 
                args=(task_queue, result_queue, self.properties)
            )
            p.start()
            processes.append(p)
            
        def producer():
            with open(self.file) as f:
                for sequence in SeqIO.parse(f, 'fasta'):
                    seq_str = str(sequence.seq).upper()
                    task_queue.put((sequence.id, seq_str))
            for _ in range(self.num_processes):
                task_queue.put(None)

        threading.Thread(target=producer).start()

        finished_processes = 0
        while finished_processes < self.num_processes:
            res = result_queue.get()
            if res is None:
                finished_processes += 1
            else:
                self.sequences_data.append(res)
                for key in self.results.keys():
                    self.results[key] += res[key]

        for p in processes:
            p.join()

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
    test_file = 'example_1000.fasta'
    
    analyzer_mp = ProteinAnalyzerMultiprocess(test_file)
    analyzer_mp.analyze()
    analyzer_mp.save_result('multiprocess_result.json')