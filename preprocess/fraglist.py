import h5py
import pandas as pd

from sys import argv

snap_file = argv[1]
out_file = argv[2]
f = h5py.File(snap_file, "r")

barcode_all = [item.decode() for item in f["BD"]['name'][:]]

barcode_UM_UQ = pd.DataFrame({'barcode': barcode_all,
                             'UM': f["BD"]['UM'][:],
                             'UQ': f["BD"]['UQ'][:]})
                             
barcode_UM_UQ.to_csv("_UM_UQ.txt", index=None, sep='\t')

frag_list = []
# bedfile = pd.DataFrame()
barcode_id = 0

for barcode in f["BD"]["name"]:  
    _chroms = f["FM"]["fragChrom"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)];
    _chroms = [item.decode() for item in _chroms];
    _start = f["FM"]["fragStart"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
    _len = f["FM"]["fragLen"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
    _barcode = [barcode.decode()] * len(_chroms)
    
    # temp = pd.DataFrame({'chr': _chroms, 
    #                  'start': _start, 
    #                  'end': _start + _len, 
    #                  'barcode': _barcode}).value_counts().rename('dup', inplace=True).reset_index()
    # bedfile = pd.concat([bedfile, temp], axis=0)
    
    frag_list_uniq = set(zip(_chroms, _start, _start + _len, _barcode)); # remove duplicated fragments
    frag_list.extend(list(frag_list_uniq))
    
    barcode_id += 1

def data_write_csv(file_name, datas):
    import codecs, csv
    file_csv = codecs.open(file_name,'w+','utf-8')
    writer = csv.writer(file_csv, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    for data in datas:
        writer.writerow(data)
    print("ok")

data_write_csv(out_file, frag_list)

# bedfile.to_csv(out_file, sep='\t', index=None)
