display_usage(){
	echo "run.sh <path to read one> <path to read two> <path to output dir>"
}
if [ $# -eq 0 ]; then
	display_usage
	exit 1
fi	
 
python3.4 pyamd.py -1 $1 -2 $2 -a ./lib/bbmap/resources/adapters.fa -r ref/mdr.fa -o $3 -b ./ref/mdr.bed -m bwa
