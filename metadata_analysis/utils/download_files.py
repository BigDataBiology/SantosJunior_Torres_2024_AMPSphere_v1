def inputsgen():
    import subprocess
    import os
    
    # create data folder
    os.makedirs('data', exist_ok=True)
    
    # publicly available in the ProGenomes resource
    print('downloading proGenomes2.1_specI_clustering.tab from ProGenomes')
    server = "progenomes.embl.de"
    address = "data"
    finput = 'proGenomes2.1_specI_clustering.tab'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}',
                     '--output-document',
                     'data/pgenomes_samples.tsv'])

    # publicly available in the GTDB resource
    print('downloading bac120_metadata_r202 from GTDB')
    server = "data.gtdb.ecogenomic.org"
    address = "releases/release202/202.0"
    finput = 'bac120_metadata_r202.tar.gz'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}',
                     '--output-document',
                     'data/bac120_metadata_r202.tar.gz'])

    # untaring file
    subprocess.call(['tar',
                     '-zxvf',
                     'data/bac120_metadata_r202.tar.gz',
                     '-C',
                     'data/'])
    
    # removing old file
    os.remove('data/bac120_metadata_r202.tar.gz')
                     
    # publicly available in the GTDB resource
    print('downloading ar122_metadata_r202 from GTDB')
    server = "data.gtdb.ecogenomic.org"
    address = "releases/release202/202.0"
    finput = 'ar122_metadata_r202.tar.gz'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}',
                     '--output-document',
                     'data/ar122_metadata_r202.tar.gz'])

    # untaring file
    subprocess.call(['tar',
                     '-zxvf',
                     'data/ar122_metadata_r202.tar.gz',
                     '-C',
                     'data/'])
    
    # removing old file
    os.remove('data/ar122_metadata_r202.tar.gz')

    # publicly available in the AMPSphere resource
    print('downloading AMPSphere_v.2021-03.fna.xz from Zenodo')
    server = "zenodo.org"
    address = "record/4606582/files"
    finput = 'AMPSphere_v.2021-03.fna.xz'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}?download=1',
                     '--output-document',
                     'data/AMPSphere_v.2021-03.fna.xz'])

    # publicly available in the AMPSphere resource
    print('downloading AMPSphere_v.2021-03.origin_samples.tsv.gz from Zenodo')
    server = "zenodo.org"
    address = "record/4606582/files"
    finput = 'AMPSphere_v.2021-03.origin_samples.tsv.gz'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}?download=1',
                     '--output-document',
                     'data/AMPSphere_v.2021-03.origin_samples.tsv.gz'])

    # generated during GMSC analysis
    print('downloading the table GMSC10.ProGenomes2.coords.txt.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/GMSC10"
    finput = "GMSC10.ProGenomes2.coords.txt.gz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading the table GMSC10.metag_smorfs.rename.txt.xz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/GMSC10"
    finput = "GMSC10.metag_smorfs.rename.txt.xz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading the table mmseqs2.lca_taxonomy.full.tsv.xz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/GMSC10"
    finput = "mmseqs2.lca_taxonomy.full.tsv.xz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during annotation of samples for AMPSphere
    print('downloading metadata')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/metadata"
    finput = 'metadata.tsv'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during Figure S1 quality analysis
    print('downloading general_envo_names.tsv')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/metadata"
    finput = 'general_envo_names.tsv'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
                     
