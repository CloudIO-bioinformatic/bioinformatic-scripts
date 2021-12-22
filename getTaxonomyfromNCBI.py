
import sys
import zipfile
import pandas as pd
from pprint import pprint
from datetime import datetime
from collections import defaultdict, Counter
from IPython.display import display

import matplotlib.pyplot as plt
plt.style.use('ggplot')

try:
    import ncbi.datasets
    assembly_accessions = []
    datasets = pd.read_csv('126_accessions.txt')
    for i in range(0, 125):
        assembly_accessions.append([datasets.values[i,j] for j in range(0, 1)])
    for data in assembly_accessions:
        print(data)
        ## start an api_instance
        api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
        genome_summary = api_instance.assembly_descriptors_by_accessions(data, limit='all')

        ## print other information
        for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
            print(
            assembly.assembly_accession,
            assembly.display_name,
            assembly.org.sci_name,
            assembly.assembly_level,
            len(assembly.chromosomes),
            assembly.submission_date,
            sep='\t')
except ImportError:
    print('ncbi.datasets module not found. To install, run `pip install ncbi-datasets-pylib`.')
