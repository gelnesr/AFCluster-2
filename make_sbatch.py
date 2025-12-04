import os
import glob
import subprocess
seed = 8

for file_ in glob.glob('/scratch/users/gelnesr/RelaxDB/fasta/*.fasta'):
    file_name = file_.split('/')[-1].strip('.fasta')
    subprocess.run(['rm', f'/scratch/users/gelnesr/tmp/afc/sbatch/{file_name}.sbatch'])

    file_new = open(f'/scratch/users/gelnesr/tmp/afc/sbatch/{file_name}.sbatch', 'x')
    file_new.write('#!/bin/bash\n')
    file_new.write(f'#SBATCH --time=2-00:00:00\n')
    file_new.write(f'#SBATCH -p possu,bioe,owners\n')
    file_new.write(f'#SBATCH --gres=gpu:1\n')
    file_new.write(f'#SBATCH --mem=40G\n')
    file_new.write(f'#SBATCH -c 8\n')
    file_new.write(f'#SBATCH --job-name={file_name}\n')
    file_new.write(f'#SBATCH --output=/scratch/users/gelnesr/tmp/afc/out/{file_name}.out\n')
    file_new.write(f'#SBATCH --error=/scratch/users/gelnesr/tmp/afc/err/{file_name}.err\n')
    file_new.write(f'#SBATCH --gpu_cmode=shared\n') 
    file_new.write(f'ml gcc/12.4.0\n')
    file_new.write(f'bash scripts/setup_slurm.sh\n')
    file_new.write(f'python afcluster.py --input {file_}')
    file_new.close()
    subprocess.run(["sbatch", f'/scratch/users/gelnesr/tmp/afc/sbatch/{file_name}.sbatch'])

