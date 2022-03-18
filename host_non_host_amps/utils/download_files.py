def inputsgen():
    import subprocess
    import os

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

    # generated during Figure 2
    print('downloading the table taxonomy_annotation.tsv')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_2/analysis/"
    finput = "taxonomy_annotation.tsv"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during annotation of samples
    print('downloading samples-min500k-assembly-prodigal-stats.tsv')
    server = "raw.githubusercontent.com"
    address = "BigDataBiology/global_data/master/freeze.v2"
    finput = 'samples-min500k-assembly-prodigal-stats.tsv'
    token = 'GHSAT0AAAAAABNWGHYYS7KTYAAA5GCIBZ2GYRTXKQQ'
    subprocess.call(['wget',
                     '-O',
                     'data/samples-min500k-assembly-prodigal-stats.tsv',
                     f'https://{server}/{address}/{finput}?token={token}'])

