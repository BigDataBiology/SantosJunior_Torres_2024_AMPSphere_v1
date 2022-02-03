# Metaproteomes quality
import gzip

# The file complete_gmsc_pgenomes_metag.tsv is available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/genes_origins
# cut -f1,2 complete_gmsc_pgenomes_metag.tsv | sed '1,1d' > ref
all_genes = set()
with open('ref', 'r') as db:
    for row in db:
        row = row.strip().split('\t')
        all_genes.add(row[0])
        
# The file metaproteome_100.tsv.gz is available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control/
green_genes = set()
with gzip.open('metaproteome_100.tsv.gz', 'rt', encoding='utf-8') as db:
    for row in db:
        if row.strip() in all_genes: green_genes.add(row.strip())

green_amps, all_amps = set(), set()
with open('ref', 'r') as db:
    for row in db:
        row = row.strip().split('\t')
        if row[0] in green_genes: green_amps.add(row[1]) 
        all_amps.add(row[1])
                
red_amps = all_amps - green_amps
with open('metaproteome_list.tsv', 'w') as out:
    out.write('AMP\tmetaproteomes\n')
    for amp in green_amps:
        out.write(f'{amp}\tgreen\n')
    for amp in red_amps:
        out.write(f'{amp}\tyellow\n')

