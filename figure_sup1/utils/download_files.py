def inputsgen():
    import subprocess
    import os

    print('generating datafolder')
    os.makedirs('data/', exist_ok=True)

    # generated during metadata analysis
    print('downloading the table gmsc_amp_genes_envohr_source.tsv.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/metadata_analysis/analysis/"
    finput = "gmsc_amp_genes_envohr_source.tsv.gz"
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
    
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
    print('downloading AMPSphere_v.2022-03.faa.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/analysis/"
    finput = 'AMPSphere_v.2022-03.faa.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading antifam results')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/data/"
    finput = 'antifam_100.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis                     
    print('downloading coordinates results')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/data/"
    finput = 'result.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading metatropeomes results')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/data/"
    finput = 'metaproteome_100.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading metatranscriptomes results')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/figure_sup1/data/"
    finput = 'metaT_100.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading negatively detected in RNAcode')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/data/"
    finput = 'rnacode_seq_false.tsv.xz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # generated during GMSC analysis
    print('downloading positively detected in RNAcode')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/data/"
    finput = 'rnacode_seq.tsv.xz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

    # publicly available in the AMPSphere resource
    print('downloading SPHERE_v.2022-03.levels_assessment.tsv.gz')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/analysis/"
    finput = 'SPHERE_v.2022-03.levels_assessment.tsv.gz'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])
 
    # generated during Figure 1 homology analysis
    print('downloading dramp annotated AMPs')
    server = "ubuntu@aws.big-data-biology.org"
    address = "/share/work/Celio/AMPSphere/v2022_03/AMPSphere_generation_v.2022-03/analysis/"
    finput = 'DRAMP_annotation.raw.tsv'
    subprocess.call(['rsync',
                     '-avzP',
                     f'{server}:{address}/{finput}',
                     'data/'])

