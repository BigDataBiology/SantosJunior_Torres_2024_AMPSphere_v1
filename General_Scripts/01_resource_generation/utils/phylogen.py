import logging
logger = logging.getLogger(__name__)


def is_command(cmds):
    """Given a command returns its path, or None.
    Given a list of commands returns the first recoverable path, or None.
    """
    try:
        from shutil import which as which  # python3 only
    except ImportError:
        from distutils.spawn import find_executable as which

    if isinstance(cmds, str):
        return which(cmds)
    else:
        for cmd in cmds:
            path = which(cmd)
            if path is not None:
                return path
        return path


def treebuilder(infile,
                ofile,
                matrix='LG+CAT',
                Nbootstraps=1000,
                stdout=None):
    """
    Builds a tree from a muscle alignment file using FastTree
    ifile - fasta file formatted alignment of peptides,
    ofile - newick tree file
    matrix - substitution matrix
    """
    import subprocess
    
    fttree_exe = is_command(['fasttree'])
    logger.debug('FastTree v2 executable: %r', fttree_exe)
    if fttree_exe is None:
        print('[ERROR] -- FastTree v2 Not Found --')
    print(f'{matrix}')
    if matrix == 'LG+CAT':
        subprocess.check_call(
            [fttree_exe,
             '-pseudo', '1',
             '-lg',
             '-cat', '20',
             '-boot', str(Nbootstraps),
             '-out', ofile,
             infile
             ], stdout=stdout)
    elif matrix == 'WAG+CAT':
        subprocess.check_call(
            [fttree_exe,
             '-pseudo', '1',
             '-wag',
             '-cat', '20',
             '-boot', str(Nbootstraps),
             '-out', ofile,
             infile
             ], stdout=stdout)
    elif matrix == 'LG':
        subprocess.check_call(
            [fttree_exe,
             '-pseudo', '1',
             '-lg',
             '-boot', str(Nbootstraps),
             '-out', ofile,
             infile
             ], stdout=stdout)
    elif matrix == 'WAG':
        subprocess.check_call(
            [fttree_exe,
             '-pseudo', '1',
             '-wag',
             '-boot', str(Nbootstraps),
             '-out', ofile,
             infile
             ], stdout=stdout)
    elif matrix == 'JTT':
        subprocess.check_call(
            [fttree_exe,
             '-pseudo', '1',
             '-fastest',
             '-out', ofile,
             infile
             ], stdout=stdout)
    else:
        print('Non-recognized matrix')


def draw_ascii(infile, ofile=None, column_width=80):
    """
    Draw an ascii-art phylogram of the given tree.
    :Parameters:

        file : file-like object
            File handle opened for writing the output drawing. (Default:
            standard output)
        column_width : int
            Total number of text columns used by the drawing.
    """
    from Bio import Phylo
    
    tree = Phylo.read(infile, 'newick')

    with open(ofile, 'w') as f:
        Phylo.draw_ascii(tree, file=f, column_width=200)

