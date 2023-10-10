# EGNC-PGO original dataset

To reproduce the result on paper, run this code: 
`cd experiments`
`../build/experiments/run-experiment -i datasets/manhattan_noisy/manhattan_noisy.irl -o ~/<some-path>/EGNC-PGO_test_environment/results/ -m risam -n 1 -d <mode>`

To plot, run this code: `./scripts/plot-traj -i datasets/manhattan_noisy/manhattan_noisy.irl -r ~/<some-path>/EGNC-PGO_test_environment/results/manhattan_noisy_risam_2023-10-10_13-50-46/ --legend`

Plotting code is from original riSAM code.