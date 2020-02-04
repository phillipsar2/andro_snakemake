import pandas as pd

list = pd.read_csv('models/bamqc_list.txt', sep = ' ', header = None, names = ['samples',$
for i in list['filename']:
    file = i + ‘/raw_data_qualimapReport/coverage_histogram.txt’
    sample = i.str.split('/').str[2].str.split('_stats').str[0] 
    data = pd.read_csv(file, sep = '\t', names = ['cov’, 'num_genomic_locations'])[1:200]
    histo_data[sample] = data['num_genomic_locations’]
histo_data = histo_data.reset_index()
histo_data.to_csv(r’reports/coverage_histogram.txt’, header = None, sep = ‘ ‘, mode = ‘a’)

