display_usage()
{
	echo "run.sh <path to experiment folder> <path to output folder>"
}
if [ $# -eq 0 ]; then
	display_usage
	exit 1
fi	
			 
python3 pyamd.py -i $1 -a ref/adapters.fa -r ref/mdr.fa -o $2 -b ref/mdr.bed -m bwa --varofint ref/Reportable_SNPs.xlsx
