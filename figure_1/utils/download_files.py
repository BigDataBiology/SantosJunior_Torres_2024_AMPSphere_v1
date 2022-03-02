def inputsgen():
    import subprocess
    import os

    print('generating datafolder')
    os.makedirs('data/', exist_ok=True)

    print('downloading the table gmsc_amp_genes_envohr_source.tsv')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/mundi_map/"
    finput = "gmsc_amp_genes_envohr_source.tsv"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
                 
    print('downloading AMPSphere_v.2021-03.faa.gz from Zenodo')
    server = "zenodo.org"
    address = "record/4606582/files"
    finput = 'AMPSphere_v.2021-03.faa.gz'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}?download=1',
                     '--output-document',
                     'data/AMPSphere_v.2021-03.faa.gz'])

    print('downloading databases for homology search')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/databases_homology"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}',
                     'data/'])

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

    print('downloading quality_candidates')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/quality_control"
    infile = ['quality_candidates.txt', 'high_quality_candidates.txt']
    for i in infile:
        subprocess.call(['rsync',
                         '-avzP',
                         f'{server}:{address}/{i}',
                         'data/'])
    
