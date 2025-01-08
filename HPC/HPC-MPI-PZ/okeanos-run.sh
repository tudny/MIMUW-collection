#!/bin/bash

cat <<EOF > matmul-run-example.batch
#!/bin/bash -l
#SBATCH --job-name matmul-run-example
#SBATCH --output "output/matmul-run-example-%j.out"
#SBATCH --error "output/matmul-run-example-%j.err"
#SBATCH --account "g96-1905"
#SBATCH --nodes 4
#SBATCH --tasks-per-node 1
#SBATCH --time 00:01:00

srun build/matmul -a data/output-mult-sq/A.csr -b data/output-mult-sq/B.csr -v -t 2D -g 100
EOF

mkdir -p output

echo "Run the following command to submit the job:"
echo "sbatch matmul-run-example.batch"
