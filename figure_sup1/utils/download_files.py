def inputsgen():
    import subprocess
    import os

    print('generating datafolder')
    os.makedirs('data/', exist_ok=True)

    # generated during metadata analysis
    print('downloading the table gmsc_amp_genes_envohr_source.tsv')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/mundi_map/"
    finput = "gmsc_amp_genes_envohr_source.tsv"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
    
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
    print('downloading AMPSphere_v.2021-03.faa.gz from Zenodo')
    server = "zenodo.org"
    address = "record/4606582/files"
    finput = 'AMPSphere_v.2021-03.faa.gz'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}?download=1',
                     '--output-document',
                     'data/AMPSphere_v.2021-03.faa.gz'])

    # generated during GMSC analysis
    print('downloading antifam results')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/quality_control"
    finput = 'antifam_100.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis                     
    print('downloading coordinates results')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/quality_control"
    finput = 'result.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading metatropeomes results')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/quality_control"
    finput = 'metaproteome_100.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading metatranscriptomes results')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/quality_control"
    finput = 'metaT_100.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading negatively detected in RNAcode')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/quality_control"
    finput = 'rnacode_seq_false.tsv.xz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading positively detected in RNAcode')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/quality_control"
    finput = 'rnacode_seq.tsv.xz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # publicly available in the AMPSphere resource
    print('downloading SPHERE_v.2021-03.levels_assessment.tsv.gz from Zenodo')
    server = "zenodo.org"
    address = "record/4606582/files"
    finput = 'SPHERE_v.2021-03.levels_assessment.tsv.gz'
    subprocess.call(['wget',
                     f'https://{server}/{address}/{finput}?download=1',
                     '--output-document',
                     'data/SPHERE_v.2021-03.levels_assessment.tsv.gz'])

    # generated during Figure 1 homology analysis
    print('downloading dramp annotated AMPs')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/homology_results"
    finput = 'dramp_candidates.txt'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during Figure 1 homology analysis
    print('downloading multi-annotated AMPs')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/files_for_figures/homology_results"
    finput = 'annotated_candidates.txt'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
                     
