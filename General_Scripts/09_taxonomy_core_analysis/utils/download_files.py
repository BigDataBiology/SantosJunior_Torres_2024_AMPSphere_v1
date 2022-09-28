def inputsgen():
    import subprocess
    import os

    print('generating datafolder')
    os.makedirs('data/', exist_ok=True)

    # publicly available in the ProGenomes resource
    print('downloading proGenomes2.1_specI_clustering.tab from ProGenomes')
    server = "progenomes.embl.de"
    address = "data"
    finput = 'proGenomes2.1_specI_clustering.tab'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}',
                     '--output-document',
                     'data/pgenomes_samples.tsv'])

    # publicly available in the AMPSphere resource
    print('downloading SPHERE_v.2022-03.levels_assessment.tsv.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/analysis/"
    finput = 'SPHERE_v.2022-03.levels_assessment.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
    
    # generated during metadata analysis
    print('downloading the table gmsc_amp_genes_envohr_source.tsv')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/metadata_analysis/outputs/"
    finput = "gmsc_amp_genes_envohr_source.tsv.gz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during Figure S1 quality analysis
    print('downloading quality families')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/analysis/"
    finput = 'quality_families.txt'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during Figure S1 quality analysis
    print('downloading quality_candidates')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/analysis/"
    finput = 'quality_candidates.txt'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during Figure S1 quality analysis
    print('downloading high_quality_candidates')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/analysis/"
    finput = 'high_quality_candidates.txt'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated manually for data curation of Prevotella
    # sequences
    print('downloading prevotella_species_list')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_2/analysis/"
    finput = 'prevotella_species_list.tsv'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated manually for data curation of Prevotella
    # sequences
    print('downloading 16S fasta for prevotella genus')
    server = "lpsn.dsmz.de"
    address = "fasta"
    finput = '37185.fasta'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}',
                     '--output-document',
                     'data/37185.fasta'])

    # generated during metadata analysis
    print('downloading the table complete_amps_associated_taxonomy.tsv.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/metadata_analysis/data/"
    finput = "complete_amps_associated_taxonomy.tsv.gz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
                     
