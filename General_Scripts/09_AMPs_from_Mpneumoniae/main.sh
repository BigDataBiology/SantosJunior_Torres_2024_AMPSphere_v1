python3 utils/download_genomes.py

mkdir analysis/

while read sample genome;
do 
	unpigz -dc data/genomes/$genome > ${genome/.gz/}
	bwa index ${genome/.gz/}
	bwa mem ${genome/.gz/} data/mgpa_ref.fa | samtools view -Sb | bedtools bamtobed > tmp.bed
	bedtools getfasta -fi ${genome/.gz/} -bed tmp.bed > analysis/${sample}.fa
	rm -rf ${genome/.gz/} *fai *pac *sa *amb *ann *bwt *bed
done < data/genomes_list

python3 utils/workseq.py

clustalo -i analysis/p1_genes.fasta \
         -o analysis/p1_genes.aln \
         -v
         
fasttree -nt \
         -gtr \
         analysis/p1_genes.aln > analysis/p1_genes.nwk

# after manual curation to get the types at iTOL, do:
#grep -f map/type1.txt ../fam_red.tsv | wc -l
#grep -f map/type2.txt ../fam_red.tsv | wc -l
#grep -f map/type3.txt ../fam_red.tsv | wc -l
