import os
import sys
import subprocess

###comment: submit information related and inputs used for running the nest program
if not os.path.exists('qsub_scripts'):
    os.mkdir('qsub_scripts')
job_ids = open('Jobid.txt', 'w')

job_nodes = 1
job_threads = 10
job_mem = '16gb'
job_walltime = '20:00:00'
job_queue = sys.argv[1]
job_email = sys.argv[2]
job_output = sys.argv[3]
file_path = os.path.dirname(os.path.abspath(__file__))

for i in range(1,102):
    qsub_file = open('qsub_scripts/qsub_nejm{0}.pbs'.format(i), 'w')
    qsub_file.write('#PBS -N Nest_nejm{0}\n#PBS -l nodes={1}:ppn={2}\n#PBS -l pmem={3}\n#PBS -l walltime={4}\n#PBS -q {5}\n#PBS -j oe\n#PBS -m abe\n#PBS -M {6}\ncd ~/\nmodule load anaconda3/latest\nsource activate nest\npython3 {7}/nest.py -i {7}/fq/NEJM/NEJM{0}.txt -r {7}/ref/mtuberculosis/mtuberculosis.fa -a {7}/ref/mtuberculosis/adapters.fa -b {7}/ref/mtuberculosis/nejm.bed -o {8}/NEJM{0} --varofint {7}/ref/mtuberculosis/nejm.csv --threads 15 --purge\n'.format(i, job_nodes, job_threads, job_mem, job_walltime, job_queue, job_email, file_path, job_output) 
    #python3 ~/projects/NeST/remover.py ~/scratch/ncbi/public/sra ~/projects/NeST/fq/NEJM/NEJM{0}.txt\n'.format(i))
    qsub_file.close()
    qrun = ['qsub', 'qsub_scripts/qsub_nejm{0}.pbs'.format(i)]
    qrunner = subprocess.Popen(qrun, stdout=subprocess.PIPE, shell=False)
    jid = qrunner.communicate()[0]
    job_ids.write('{0}'.format(jid.decode('UTF-8')))

job_ids.close()
