def make_seqlogo():
    '''
    Generate sequence logo per alignment in the 
    selected families
    
    :input: aligned peptides fasta file
    
    :output: pdf with the seqlogo
    
    Dependencies = Weblogo
    
    https://github.com/WebLogo/weblogo/
    '''
    import os
    import glob
    import subprocess
    
    data_folder = 'analysis/families/aln'
    analysis_folder = 'analysis/families/seqlogo'
    
    os.makedirs(analysis_folder, exist_ok=True)
    
    for f in glob.glob(f'{data_folder}/*.aln'):

        fasta_name = f.split('/')[-1].replace('.aln', '')

        seqlogo_out = f"{analysis_folder}/{fasta_name}.pdf"

        subprocess.call(['weblogo',
                         '-f', f,
                         '-D', 'fasta',
                         '-A', 'protein',
                         '-o', seqlogo_out,
                         '-F', 'pdf',
                         '-U', 'bits',
                         '--composition', 'auto',
                         '-t', fasta_name,
                         '-n', '100',
                         '--number-interval', '1', 
                         '--errorbars', 'NO',
                         '--fontsize', '36',
                         '--number-fontsize', '24',
                         '--stack-width', '50',
                         '--box', 'NO',
                         '--resolution', '300',
                         '--scale-width', 'NO'])
                         
        print(seqlogo_out)
       

