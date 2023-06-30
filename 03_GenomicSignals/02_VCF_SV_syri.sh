cd /home/jnrunge/data/trd/Victor_SV/alignments

strain1=$1
strain2=$2

echo -e "/home/vloegler/SVCalling2/03-data/BestCollapsed/$strain1.Final.fasta\t$strain1" > "${strain1}_${strain2}.mi80l100syri.genomes"
echo -e "/home/vloegler/SVCalling2/03-data/BestCollapsed/$strain2.Final.fasta\t$strain2" >> "${strain1}_${strain2}.mi80l100syri.genomes"

grep "^>" /home/vloegler/SVCalling2/03-data/BestCollapsed/$strain1.Final.fasta | cut -f 1 -d' ' | cut -f 2 -d'>' > "${strain1}_${strain2}.mi80l100syri.chrs"

nucmer -p ${strain1}_${strain2} /home/vloegler/SVCalling2/03-data/BestCollapsed/$strain1.Final.fasta /home/vloegler/SVCalling2/03-data/BestCollapsed/$strain2.Final.fasta

delta-filter -m -i 80 -l 100 ${strain1}_${strain2}.delta > ${strain1}_${strain2}.mi80l100.delta

show-coords -THrd ${strain1}_${strain2}.mi80l100.delta > ${strain1}_${strain2}.mi80l100.delta.coords

syri -c ${strain1}_${strain2}.mi80l100.delta.coords -d ${strain1}_${strain2}.mi80l100.delta -r /home/vloegler/SVCalling2/03-data/BestCollapsed/$strain1.Final.fasta -q /home/vloegler/SVCalling2/03-data/BestCollapsed/$strain2.Final.fasta --prefix ${strain1}_${strain2}.mi80l100

plotsr --sr ${strain1}_${strain2}.mi80l100syri.out -o ${strain1}_${strain2}.mi80l100syri.out.itx.pdf --genomes ${strain1}_${strain2}.mi80l100syri.genomes --itx --chrord ${strain1}_${strain2}.mi80l100syri.chrs -H 5 -W 5

plotsr --sr ${strain1}_${strain2}.mi80l100syri.out -o ${strain1}_${strain2}.mi80l100syri.out.pdf --genomes ${strain1}_${strain2}.mi80l100syri.genomes --chrord ${strain1}_${strain2}.mi80l100syri.chrs

gzip -f ${strain1}_${strain2}.mi80l100syri.out

. ~/activate.sh bwaetc

bgzip -f ${strain1}_${strain2}.mi80l100syri.vcf

bcftools index ${strain1}_${strain2}.mi80l100syri.vcf.gz

