srun -w c1 -n32 -p bcc --pty bash
cd /net/bmc-lab3/data/bcc/projects/cmelen-Love/ont_flye_workflow_runs

git clone -b sequential_v3 https://github.com/KochInstitute-Bioinformatics/ont-flye-workflow.git
cd ont-flye-workflow
git branch

git checkout -b sequential_v4
git push -u origin sequential_v4

module add miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
#conda activate nf-core_v24
conda activate nf-core_v25
module add singularity/3.10.4
