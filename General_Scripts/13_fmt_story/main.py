from utils.fmt_startup import fmt_startup
from utils.evaluation import eval_results
from utils.merger import merger
from utils.ftests import ftests
from utils.plot_per_genome import plot_genomes

def main():
    fmt_startup()
    eval_results()
    merger()
    ftests()
    plot_genomes()
        
    
if __name__ == '__main__':
    main()
    
