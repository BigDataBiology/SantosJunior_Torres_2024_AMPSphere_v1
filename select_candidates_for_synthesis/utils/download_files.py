def inputsgen():
    import subprocess
    import os
    
    # create data folder
    os.makedirs('data', exist_ok=True)
    
    # publicly available in the AMPSphere resource
    print('downloading SPHERE_v.2022-03.levels_assessment.tsv.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/analysis/"
    finput = 'SPHERE_v.2022-03.levels_assessment.tsv.gz'
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

    # internal resource in AMPSphere
    print('downloading high_quality_candidates.txt')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/analysis/"
    finput = 'high_quality_candidates.txt'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # internal resource in AMPSphere
    print('downloading quality_candidates.txt')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/analysis/"
    finput = 'quality_candidates.txt'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

