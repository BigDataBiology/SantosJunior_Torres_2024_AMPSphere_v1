import gzip

print('catching green genes')
# The file result.tsv.gz is available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control/
# zgrep "T" result.tsv.gz | pigz --best > true_result.tsv.gz
green_genes = set()
with open('true_genes_in_AMPSphere.txt', 'r') as db:
    for idx, row in enumerate(db):
        green_genes.add(row.strip().split('\t')[0])
        print(idx)

print('getting amps')
# The file complete_gmsc_pgenomes_metag.tsv is available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/genes_origins
# cut -f1,2 complete_gmsc_pgenomes_metag.tsv | sed '1,1d' > ref
green_amps, all_amps = set(), set()
with open('ref', 'r') as db:
    for row in db:
        row = row.strip().split('\t')
        if row[0] in green_genes: green_amps.add(row[1]) 
        all_amps.add(row[1])

print('getting reds')                
red_amps = all_amps - green_amps
with open('coordinates_list.tsv', 'w') as out:
    out.write('AMP\tcoordinates\n')
    for amp in green_amps:
        out.write(f'{amp}\tgreen\n')
    for amp in red_amps:
        out.write(f'{amp}\tred\n')

