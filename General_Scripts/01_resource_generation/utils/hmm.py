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


def hbuild(ifile, ofile, fam, stdout=None):
    """
    Function to perform create HMM profiles. It uses the HMMer
    hmm build software and the alignment of each family of
    peptides.
    """
    import subprocess
    
    hmmbuild_exe = is_command(['hmmbuild'])
    logger.debug('HMMbuild executable: %r', hmmbuild_exe)
    if hmmbuild_exe is None:
        print('[ERROR] -- HMMbuild Not Found --')
    subprocess.check_call(
        [hmmbuild_exe,
         '--amino',
         '-n', fam,
         ofile,
         ifile],
        stdout=stdout)


def clean_df(df):
    """
    Remove nan columns and remaining rows containing nan values.
    """
    df = df.dropna(axis='columns', how='all')
    df = df.dropna(axis='rows', how='any')
    df = df.astype('float64')
    df.index = df.index.astype('int64')
    return df


def get_indelinfo(hmmlogo_df):
    """
    Return indelinfo for hidden Markov model logo.
    """
    selection = hmmlogo_df['Y'].isna()
    indelinfo = hmmlogo_df[selection].copy()
    indelinfo = clean_df(indelinfo)
    indelinfo = indelinfo.rename(columns=dict(A='insert probability',
                                              C='average insert length',
                                              D='occupancy'))
    return indelinfo


def read_hmmlogo(hmmlogo_output):
    """
    Return heights and indelinfo from hidden Markov model logo.
    """
    from io import BytesIO
    from pandas import read_csv
    
    column_names = ['profile']+list('ACDEFGHIKLMNPQRSTVWY-')
    seperator = ':?\s+\(?\s?'
    hmmlogo_df = read_csv(
                          BytesIO(hmmlogo_output),
                          index_col=0,
                          skiprows=2,
                          header=None,
                          sep=seperator,
                          names=column_names,
                          engine="python",
                          )
    hmmlogo_df = hmmlogo_df.drop('-', axis='columns')
    heights = clean_df(hmmlogo_df)
    indelinfo = get_indelinfo(hmmlogo_df)
    return heights, indelinfo


def get_hmmlogo(hmm):
    """
    Return hidden Markov model logo from hidden Markov model.
    """
    from subprocess import run
    
    in_file = open(hmm, "rb")
    data = in_file.read()
    in_file.close()

    hmmlogo_exe = is_command(['hmmlogo'])
    logger.debug('HMMlogo executable: %r', hmmlogo_exe)
    if hmmlogo_exe is None:
        print('[ERROR] -- HMMlogo Not Found --')

    arguments = [hmmlogo_exe, '/dev/stdin']
    process = run(arguments, input=data, capture_output=True)
    heights, indelinfo = read_hmmlogo(process.stdout)
    return heights, indelinfo


def initialize_plot(heights):
    import matplotlib.pyplot as plt

    length = len(heights)
    fig, ax = plt.subplots(figsize=(length/2, 7))
    ax.set_xlim(0, length)
    return fig, ax


def plot_heights(alph, cmap, heights, ax):
    import matplotlib.pyplot as plt

    y_max = 0
    for x_offset, (_, profile_heights) in enumerate(heights.iterrows()):
        plot_profile(alph, cmap, profile_heights, x_offset, ax)
        y_max = max(y_max, profile_heights.sum())
    ax.set_ylim(0, y_max)


def get_aa_patch(alph, cmap, aminoacid, scale, offset):
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    
    path = Path(*alph[aminoacid])
    path.vertices *= scale
    path.vertices += offset
    patch = PathPatch(path, color=cmap[aminoacid])
    return patch


def plot_profile(alph, cmap, profile, x_offset, ax):
    y_offset = 0
    profile = profile.sort_values()
    for aminoacid, height in profile.iteritems():
        patch = get_aa_patch(alph,
                             cmap,
                             aminoacid,
                             [1, height],
                             [x_offset, y_offset])
        ax.add_patch(patch)
        y_offset += height


def get_svg(fig, ax, indelinfo, ofile):
    from io import StringIO
    import matplotlib.pyplot as plt

    ax.xaxis.set_visible(False)
    indelinfo = indelinfo.round(decimals=2).T
    ax.table(
             cellText=indelinfo.values,
             rowLabels=indelinfo.index,
             colLabels=indelinfo.columns,
             )
    fig.tight_layout()
    fig.savefig(ofile, format='svg', bbox_inches='tight')
    plt.close()
    

def pathdicts(alph, cmap):
    '''
    Loads the color files
    '''
    import json
    
    with open(alph, 'r') as f:
        alph = json.load(f)
    with open(cmap, 'r') as f:
        cmap = json.load(f)

    return alph, cmap

def pict_hmmlogo(alph, cmap, ifile, ofile):
    alph, cmap = pathdicts(alph, cmap)
    heights, indelinfo = get_hmmlogo(ifile)
    fig, ax = initialize_plot(heights)
    plot_heights(alph, cmap, heights, ax)
    get_svg(fig, ax, indelinfo, ofile)

