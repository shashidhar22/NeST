time ./pyamd.py -i fq/test_dataset -r ref/mdr.fa -a ref/adapters.fa -b ref/mdr.bed -o local/experiment1/ -m bwa --varofint ref/Reportable_SNPs_Report_v2.xlsx
echo "# Experiment1\nRun MaRS on dataset13 using BWA aligner" > local/experiment1/details.md
#time ./pyamd.py -i fq/test_dataset -r ref/mdr.fa -a ref/adapters.fa -b ref/mdr.bed -o local/experiment2/ -m bowtie2 --varofint ref/Reportable_SNPs_Report_v2.xlsx
#echo "# Experiment2\nRun MaRS on dataset13 using Bowtie2 aligner" > local/experiment2/details.md
#time ./pyamd.py -i fq/test_dataset -r ref/mdr.fa -a ref/adapters.fa -b ref/mdr.bed -o local/experiment3/ -m bbmap --varofint ref/Reportable_SNPs_Report_v2.xlsx
#echo "# Experiment3\nRun MaRS on dataset13 using BBMap aligner" > local/experiment3/details.md
#time ./pyamd.py -i fq/test_dataset -r ref/mdr.fa -a ref/adapters.fa -b ref/mdr.bed -o local/experiment4/ -m snap --varofint ref/Reportable_SNPs_Report_v2.xlsx
#echo "# Experiment4\nRun MaRS on dataset13 using SNAP aligner" > local/experiment4/details.md

