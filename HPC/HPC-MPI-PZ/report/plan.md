# Matmul

## Testing plan

Okeanos setup enables 24 cores per node, and let's assume 8 nodes at max per run.
Matrix have to be squared.

$$ N \in \{ 1000, 10000, 10000 \} $$ 
$$ d \in \{ 0.3,  0.01,  0.05  \} $$

$$ p \in \{ (24 * 6, 1), (24 * 6, 4), (24 * 6, 9) \} $$

$$ alg \in \{ 2D, 3D \} $$
