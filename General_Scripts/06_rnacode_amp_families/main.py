from utils.genes_to_clusters import merge_clusters, create_fasta
from utils.rnacode_call import batch_process
from utils.largefams import large_fam_clusters, fasta_fragments
from utils.rnacode_lfam import frag_rnacode
from utils.lfam_arrange import lfam_rnacode_output

def main():
    create_fasta(merge_clusters())
    batch_process()
    fasta_fragments(large_fam_clusters())
    frag_rnacode()
    lfam_rnacode_output()
    

if __name__=='__main__':
    main()

