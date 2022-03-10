from __future__ import with_statement

import base64
import itertools
import logging
import os
import subprocess


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


def hbuild(ifile, ofile, fam, stdout=None):
    """
    Function to perform create HMM profiles. It uses the HMMer
    hmm build software and the alignment of each family of
    peptides.
    """
    hmmbuild_exe = is_command(['hmmbuild'])
    logger.debug('HMMbuild executable: %r', hmmbuild_exe)
    if hmmbuild_exe is None:
        print('[ERROR] -- HMMbuild Not Found --')
    subprocess.check_call(
        [hmmbuild_exe,
         '-n', fam,
         ofile,
         ifile],
        stdout=stdout)


def call_mmseqs(fin, fout, db, threads, tmp='tmp', stdout=None):
    """
    Function to perform peptides searching using MMSeqs-2. Fixed
    parameters include: global identity of 90%, query coverage
    of 50%, output in m8 format and usage of a tmp temporary
    folder.
    """
    mmseq_exe = is_command(['mmseqs'])
    logger.debug('MMSeqs executable: %r', mmseq_exe)
    if mmseq_exe is None:
        print('[ERROR] -- MMSeqs Not Found --')
    subprocess.check_call(
        [mmseq_exe,
         'easy-search',
         fin,
         db,
         fout,
         tmp,
         '--threads', str(threads)],
        stdout=stdout)


def call_cdhit(fin, fout, threshold, wordsize, threads, stdout=None):
    """
    Function to perform peptides clustering using cd-hit. Fixed
    parameters include: global identity, print alignment overlap,
    cluster sorting by decreasing size in the clstr and fasta
    files, slow clustering mode, band_width (5), throw away
    sequences under 6 -- safe margin, minimum coverage of 80% of
    the shorter sequence, print entire access code until first
    space, maximum disk memory usage of 16GB.
    """
    cdhit_exe = is_command(['cd-hit', 'cdhit'])
    logger.debug('Cdhit executable: %r', cdhit_exe)
    if cdhit_exe is None:
        print('[ERROR] -- Cdhit Not Found --')
    subprocess.check_call(
        [cdhit_exe,
         '-i', fin,
         '-o', fout,
         '-c', str(threshold),
         '-n', str(wordsize),
         '-G', '1',
         '-g', '1',
         '-b', '5',
         '-l', '6',
         '-p', '1',
         '-sf', '1',
         '-sc', '1',
         '-aS', '0.8',
         '-M', '16000',
         '-d', '0',
         '-T', str(threads)],
        stdout=stdout)


def call_crev(fin, fout, stdout=None):
    """
    Function to combine a .clstr file with its parent .clstr file
    after a cd-hit hierarchical clustering procedure. It uses the
    script clstr_rev.pl in the PATH as clstr_rev.
    """
    crev_exe = is_command(['clstr_rev'])
    if crev_exe is None:
        print('[ERROR] -- Clstr_rev Not Found --')
    subprocess.check_call([crev_exe, fin, fout], stdout=stdout)


def delimit(iterable, splitstring):
    """
    Function to split cd-hit cluster files in a sensible way.
    Groups the items in the list using a keyword, and creates a
    list of sublists that are grouped accordingly. Returns the
    list of sublists.
    """
    return [list(g) for k, g in itertools.groupby(
        iterable, lambda x:x in splitstring) if not k]

