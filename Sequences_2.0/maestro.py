import os  
from multiprocessing import Pool  
processes = ('main1.py', 'main2.py', 'main3.py')  
def run_process(process):  
 os.system('python {}'.format(process))  
pool = Pool(processes=3)  
pool.map(run_process, processes) 