def file_size(infile):
    import os
    file_stats = os.stat(infile)
    return file_stats.st_size / (1024 * 1024)


def open_file(infile):
    if infile.endswith('.tar.xz'):
        import tarfile
        tarf = tarfile.open(infile, 'r:xz')
        x = tarf.getnames()
    elif infile.endswith('gz'):
        import gzip
        x = gzip.open(infile, 'rt')
        x = x.readlines()[:5]
    elif infile.endswith('xz'):
        import lzma
        x = lzma.open(infile, 'rt')
        x = x.readlines()[:5]
    else:
        return 'file format not known'
    return x
   
        
def hashfile(infile):
    import hashlib as hash
    BLOCKSIZE=65536
    md5 = hash.md5()
    with open(infile, 'rb') as ifile:
        file_buffer = ifile.read(BLOCKSIZE)
        while len(file_buffer) > 0:
            md5.update(file_buffer)
            file_buffer = ifile.read(BLOCKSIZE)
    return md5.hexdigest()
   
    
def writeoutput(infile):
    print(f'Processing file {infile}')
    md5 = hashfile(infile)
    header = open_file(infile)
    size = file_size(infile)
    print('Outputting results')
    title = '.'.join(infile.split('.')[:-1])
    name = infile.split('/')[-1]
    with open(f'docs/{name}.md', 'wt') as ofile:
        ofile.write(f'**{title}**\n\n')
        ofile.write(f'**Description:**\n\n')
        ofile.write(f'**MD5 SUM:**\t{md5}\n\n')
        ofile.write(f'**Size (MBytes):**\t{size}\n\n')
        ofile.write('**Content sample (first 5 items):**\n\n')
        for i in header:
            if i.endswith('\n'): ofile.write(i)
            else: ofile.write(f'{i}\n')
        ofile.write('[...]\n')
        

def main():
    from glob import glob
    for infile in glob('*z'):
        writeoutput(infile)


if __name__ == '__main__':
    main()
    
