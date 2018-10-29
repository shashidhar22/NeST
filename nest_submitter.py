import os
import sys
import subprocess

if not os.path.exists('qsub_scripts'):
    os.mkdir('qsub_scripts')
job_ids = open('Jobid.txt', 'w')

for i in range(1,102):
    qsub_file = open('qsub_scripts/qsub_nejm{0}.pbs'.format(i), 'w')
    qsub_file.write('#PBS -N Nest_nejm{0}\n#PBS -l nodes=1:ppn=10\n#PBS -l pmem=16gb\n#PBS -l walltime=20:00:00\n#PBS -q biocluster-6\n#PBS -j oe\n#PBS -m abe\n#PBS -M pace.shashidhar@gmail.com\ncd ~/\nmodule load anaconda3/latest\nsource activate nest\npython3 ~/projects/NeST/nest.py -i ~/projects/NeST/fq/NEJM/NEJM{0}.txt -r ~/projects/NeST/ref/mtuberculosis/mtuberculosis.fa -a ~/projects/NeST/ref/mtuberculosis/adapters.fa -b ~/projects/NeST/ref/mtuberculosis/nejm.bed -o ~/scratch/Nest/results/NEJM{0} --varofint ~/projects/NeST/ref/mtuberculosis/nejm.csv --threads 15 --purge\npython3 ~/projects/NeST/remover.py ~/scratch/ncbi/public/sra ~/projects/NeST/fq/NEJM/NEJM{0}.txt\n'.format(i))
    qsub_file.close()
    qrun = ['qsub', 'qsub_scripts/qsub_nejm{0}.pbs'.format(i)]
    qrunner = subprocess.Popen(qrun, stdout=subprocess.PIPE, shell=False)
    jid = qrunner.communicate()[0]
    job_ids.write('{0}'.format(jid.decode('UTF-8')))

job_ids.close()
