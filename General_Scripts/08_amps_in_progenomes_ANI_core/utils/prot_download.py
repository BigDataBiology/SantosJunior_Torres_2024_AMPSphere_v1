def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    from Bio import Entrez
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle, validate=False)
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
    Entrez.email = "A.N.Other@example.com"
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
        link = os.path.join(url,label+'_protein.faa.gz')
        print (link)
        links.append(link)
        if download == True:
            #download link
            urllib.request.urlretrieve(link, f'{label}.faa.gz')
            os.rename(f'{label}.faa.gz', f'{term}.faa.gz')
    return links


def batchdownload(samplelist):
    for x in samplelist:
        print(x)
        try:
            links = get_assemblies(x,
                                   download=True)
        except:
            print(f'Wrong link to download -- ERROR')
        
