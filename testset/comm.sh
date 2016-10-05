kmer=31
meta_p=../bin/metacortex_k${kmer}

export R_ENV_PATH='../'

b=100
n=15

cortex_file='all.ctx'
contig_file='contigs.fa'
log_file='log.txt'
file_list='allfiles.txt'

if [ -f ${cortex_file} ] ; then
	rm ${cortex_file}
fi
if [ -f ${contig_file} ] ; then
        rm ${contig_file}
fi

mkdir graphs; chmod 755 graphs

echo 'basic_genome.fa' > ${file_list}

${meta_p} -k ${kmer} -n ${n} -b ${b} -i ${file_list} -t fasta -o ${cortex_file} -f ${contig_file} -g 100 -l ${log_file}  -S
