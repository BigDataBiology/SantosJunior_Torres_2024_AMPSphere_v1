import gzip

# The file antifam_100.tsv.gz available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control/
red_genes = set()
with gzip.open('antifam_100.tsv.gz', 'rt', encoding='utf-8') as db:
    for row in db:
        red_genes.add(row.strip())

# The file complete_gmsc_pgenomes_metag.tsv is available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/genes_origins
# cut -f1,2 complete_gmsc_pgenomes_metag.tsv | sed '1,1d' > ref
red_amps, all_amps = set(), set()
with open('ref', 'r') as db:
    for row in db:
        row = row.strip().split('\t')
        if row[0] in red_genes: red_amps.add(row[1]) 
        all_amps.add(row[1])
                
green_amps = all_amps - red_amps

with open('antifam_list.tsv', 'w') as out:
    out.write('AMP\tAntifam\n')
    for amp in green_amps:
        out.write(f'{amp}\tgreen\n')
    for amp in red_amps:
        out.write(f'{amp}\tred\n')

