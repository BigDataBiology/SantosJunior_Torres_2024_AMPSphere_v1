def inputsgen():
    import subprocess
    import os
    import shutil

    print('generating datafolder')
    os.makedirs('data/', exist_ok=True)

    # generated during metadata analysis
    print('downloading the table gmsc_amp_genes_envohr_source.tsv')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/metadata_analysis/outputs/"
    finput = "gmsc_amp_genes_envohr_source.tsv.gz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
                 
    # publicly available in the AMPSphere resource
    print('downloading AMPSphere_v.2022-03.faa.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/analysis/"
    finput = 'AMPSphere_v.2022-03.faa.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # publicly available in different databases
    print('downloading databases for homology search')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_1/data/databases_homology"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}',
                     'data/'])
    
    # moving files from inside folders
    # substituting internal soflinks
    truepep = 'true_pep_2022_vs_progenomesgmgc.tsv.xz'
    conver = 'data/databases_homology/converter'
    os.remove(f'data/databases_homology/{truepep}')
    os.rename(f'{conver}/analysis/{truepep}',
              f'data/databases_homology/{truepep}')              
    shutil.rmtree(conver)
    
    # GMGC resource
    print('downloading GMGC')
    ask = input('Do you want to download the GMGC catalog? It may consume much space! (y/n)')
    if (ask == 'Y') or (ask == 'y'):
        server = "ubuntu@aws.big-data-biology.org"
        address = "/share/resources/GMGC10"
        finput = "GMGC10.proGenomes.faa"
        subprocess.call(['rsync',
                         '-avzP',
                         f'{server}:{address}/{finput}',
                         'data/databases_homology/'])
    
    # from the AMPSphere resource, generated during computations for Sup. Fig. 1
    print('downloading quality_candidates')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/analysis/"
    infile = ['quality_assessment.tsv', 'quality_candidates.txt', 'high_quality_candidates.txt']
    for i in infile:
        subprocess.call(['rsync',
                         '-avzP',
                         f'{server}:{address}/{i}',
                         'data/'])
    
