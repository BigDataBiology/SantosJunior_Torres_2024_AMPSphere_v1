import lzma

print('generating all genes')
all_genes = set()
with open('ref', 'r') as db:
    for row in db:
        row = row.strip().split('\t')
        all_genes.add(row[0])

print('generating red genes')
red_genes = set()
with lzma.open('rnacode_seq_false.tsv.xz', 'rt', encoding='utf-8') as db:
    for row in db:
        if row.strip() in all_genes: red_genes.add(row.strip())

# The file rnacode_seq.tsv.xz available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control/
print('generating green genes')
green_genes = set()
with lzma.open('rnacode_seq.tsv.xz', 'rt', encoding='utf-8') as db:
    for row in db:
        if row.strip() in all_genes: green_genes.add(row.strip())

# The file complete_gmsc_pgenomes_metag.tsv is available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/genes_origins
# cut -f1,2 complete_gmsc_pgenomes_metag.tsv | sed '1,1d' > ref
print('generating AMP lists')
red_amps, green_amps, all_amps = set(), set(), set()
with open('ref', 'r') as db:
    for row in db:
        row = row.strip().split('\t')
        all_amps.add(row[1])
        if row[0] in red_genes: red_amps.add(row[1])
        elif row[0] in green_genes: green_amps.add(row[1])

print('generating yellow AMPs')
yellow_amps = all_amps - green_amps 
yellow_amps = yellow_amps - red_amps

print('writing out')                
with open('rnacode_list.tsv', 'w') as out:
    out.write('AMP\tRNAcode\n')
    for amp in green_amps:
        out.write(f'{amp}\tgreen\n')
    for amp in red_amps:
        out.write(f'{amp}\tred\n')
    for amp in yellow_amps:
        out.write(f'{amp}\tyellow\n')

