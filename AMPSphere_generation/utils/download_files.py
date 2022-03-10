def inputsgen():
    import subprocess
    import os
    
    # create data folder
    os.makedirs('data', exist_ok=True)
    
    # patch colors for the HMMlogo package
    print('downloading ALPHA')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/data"
    finput = 'alph.json'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # patch colors for the HMMlogo package
    print('downloading CMAP')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/data"
    finput = 'cmap.json'
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
                     'data/pgenomes_samples.tsv'])

    # publicly available in the ProGenomes resource
    print('downloading proGenomes2.1_specI_lineageNCBI.tab from ProGenomes')
    server = "progenomes.embl.de"
    address = "data"
    finput = 'proGenomes2.1_specI_lineageNCBI.tab'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}',
                     '--output-document',
                     'data/proGenomes2.1_specI_lineageNCBI.tab'])

    # resource from ProGenomes v2 with the genomes
    # filtered in detail withtout those black listed
    print('downloading proGenomes2.1 blacklist')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/data"
    finput = 'combinedVsearch3.50p.20n.norm.opt.0.042.filteredBlackList.map'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    # AMP predictions for metaG and genomes
    # (GMSC < 34M indicate prediction from Genomes)
    print('downloading the table GMSC10.Macrel_05.AMPs.tsv.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/GMSC10"
    finput = "GMSC10.Macrel_05.AMPs.tsv.gz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    # contains the gene features for all metaG derived 
    # smorfs
    print('downloading the table GMSC10.metag_smorfs.rename.txt.xz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/GMSC10"
    finput = "GMSC10.metag_smorfs.rename.txt.xz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # resource available publicly in DRAMP v2 website
    print('downloading databases for homology search')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/databases_homology/"
    finput = "DRAMP.fa"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # resource manually curated from the data
    # available publicly in DRAMP v2 website
    print('downloading databases for homology search')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/data"
    finput = "dramp_anno.py"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # pre-computed resource generated during AMPSphere v.2021-03
    # it was generated using the original mmseqs version
    # because of that results are discrepant from recent 
    # versions
    print('downloading the table DRAMP_filter.raw.tsv')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/data"
    finput = "DRAMP_filter.raw.tsv"
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
    print('downloading the table GMSC10.proGenomesv2_05.AMPs.rename.tsv.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/GMSC10"
    finput = "GMSC10.proGenomesv2_05.AMPs.rename.tsv.gz"
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
    # genes selected from the large file
    # generated elsewhere
    print('downloading the table gmsc_genes.fna.xz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/data"
    finput = "gmsc_genes.fna.xz"
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

