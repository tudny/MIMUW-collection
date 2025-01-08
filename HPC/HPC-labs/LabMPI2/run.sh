#!/bin/bash

SIZES=(100 1000 10000 100000)

for nodes in {2..16}; do
    for size in "${SIZES[@]}"; do
cat <<EOF >floyd-warshall-"$nodes"-"$size".batch
#!/bin/bash -l
#SBATCH --job-name floyd-warshall-$nodes-$size
#SBATCH --output "floyd-warshall-$nodes-$size-%j.out"
#SBATCH --error "floyd-warshall-$nodes-$size-%j.err"
#SBATCH --account "g96-1905"
#SBATCH --nodes $nodes
#SBATCH --tasks-per-node 1
#SBATCH --time 00:05:00

srun floyd-warshall-par.exe --show-results $size
EOF
    echo floyd-warshall-"$nodes"-"$size".batch
    done
done
