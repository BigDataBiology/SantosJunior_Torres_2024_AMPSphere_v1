def inputsgen():
    import subprocess
    import os
    
    # create data folder
    os.makedirs('data', exist_ok=True)
    
    # publicly available in the ProGenomes resource
    print('downloading proGenomes2.1_specI_lineageNCBI.tab from ProGenomes')
    server = "progenomes.embl.de"
    address = "data"
    finput = 'proGenomes2.1_specI_lineageNCBI.tab'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}',
                     '--output-document',
                     'data/proGenomes2.1_specI_lineageNCBI.tab'])

    # internal resource // black listed taxonomies
    print('downloading lost lineages')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/metadata_analysis/data/"
    finput = 'ncbi_missing_lineages.txt'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
    
    # publicly available in the ProGenomes resource
    print('downloading proGenomes2.1_specI_clustering.tab from ProGenomes')
    server = "progenomes.embl.de"
    address = "data"
    finput = 'proGenomes2.1_specI_clustering.tab'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}',
                     '--output-document',
                     'data/progenomes_samples.tsv'])

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
    print('downloading AMPSphere_v.2022-03.fna.xz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/analysis/"
    finput = 'AMPSphere_v.2022-03.fna.xz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # publicly available in the AMPSphere resource
    print('downloading AMPSphere_v.2022-03.origin_samples.tsv.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/analysis/"
    finput = 'AMPSphere_v.2022-03.origin_samples.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

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
    server = "raw.githubusercontent.com"
    address = "BigDataBiology/global_data/master/freeze.v2"
    finput = 'metadata.tsv'
    token = 'GHSAT0AAAAAABNWGHYZ4VSZJ2CRBCOPV5BGYRPXCSA'
    subprocess.call(['wget',
                     '-O',
                     'data/metadata.tsv',
                     f'https://{server}/{address}/{finput}?token={token}'])

    # generated during annotation of samples for AMPSphere
    print('downloading general envo names')
    server = "raw.githubusercontent.com"
    address = "BigDataBiology/global_data/master/freeze.v2"
    finput = 'general_envo_names.tsv'
    token = 'GHSAT0AAAAAABNWGHYZV7HGOUYSQX35NW2KYRPX3CA'
    subprocess.call(['wget',
                     '-O',
                     'data/general_envo_names.tsv',
                     f'https://{server}/{address}/{finput}?token={token}'])

