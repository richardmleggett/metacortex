kmer=31
meta_p=../bin/metacortex_k${kmer}

export R_ENV_PATH='../'

b=100
n=15

#filename="reads/basic_genome.fa"
filename="reads/single_snp_genome.fa"
#filename="reads/single_insert_genome.fa"
#filename="reads/single_del_genome.fa"

name=`basename ${filename}`
name=`echo ${name%.*}`

cortex_file="${name}.ctx"
contig_file="${name}.fa"
log_file="${name}.txt"
file_list="allfiles.txt"

if [ -f ${cortex_file} ] ; then
	rm ${cortex_file}
fi
if [ -f ${contig_file} ] ; then
        rm ${contig_file}
fi

mkdir graphs; chmod 755 graphs

echo ${filename} > ${file_list}

${meta_p} -k ${kmer} -n ${n} -b ${b} -i ${file_list} -t fasta -o ${cortex_file} -f ${contig_file} -g 100 -l ${log_file}  -S -G graph_out.gv

dot -Tdot graph_out.gv -o graph_out.dot


