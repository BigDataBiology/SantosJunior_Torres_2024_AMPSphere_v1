def download_genomes():    
    '''
    Download genomes from the specI cluster 1969
    in ProGenomes v2 database. It uses a NCBI download
    engine that relies on the user e-mail,
    please check out through the code where to add it.
    '''
    # get samples
    import pandas as pd
    samples = pd.read_table('data/progenomes_cluster1969.tsv', sep='\t', header='infer')
    
    def get_assembly_summary(id):
        """Get esummary for an entrez id"""
        from Bio import Entrez
        esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
        esummary_record = Entrez.read(esummary_handle)
        return esummary_record
    
    
    def get_assemblies(term, download=True, path='assemblies'):
        """Download genbank assemblies for a given search term.
        Args:
            term: search term, usually organism name
            download: whether to download the results
            path: folder to save to
        """
        import os, urllib
        from Bio import Entrez
        #provide your own mail here
        Entrez.email = "A.N.Other@example.com"  # <<<<- ADD HERE USER's EMAIL
        handle = Entrez.esearch(db="assembly", term=term, retmax='200')
        record = Entrez.read(handle)
        ids = record['IdList']
        print (f'found {len(ids)} ids')
        links = []
        for id in ids:
            #get summary
            summary = get_assembly_summary(id)
            #get ftp link
            url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
            if url == '':
                continue
            label = os.path.basename(url)
            #get the fasta link - change this to get other formats
            link = os.path.join(url,label+'_genomic.fna.gz')
            print (link)
            links.append(link)
            if download == True:
                #download link
                urllib.request.urlretrieve(link, f'{label}.fna.gz')
        return links
    
    
    for x in samples['sample']:
        print(x)
        links = get_assemblies(x,
                               download=True)
    
    
