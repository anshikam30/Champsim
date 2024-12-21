##ChampSim

This is part of Assignment 1 Q2 of course E0-243 Computer architecture offered at CSA department at Indian Institute of Science.

In this assignment we evaluate different branch predictors with ChampSim simulator.

We have only modified the files in branch folder (added Tage and Hybrid Predictors).
ChampSim installation guilde can be found in ChampSim's README.md file. Follow the instructions on ChampSim Github link and set it up and the replace the Champsim's branch folder with the given branch folder, and then run.



## Process:

We have evaluated three cloudsuite traces  `cloud9_phase0_core1.trace.xz, nutch_phase1_core0.trace.xz , streaming_phase1_core0.trace.xz`,and one SPEC trace `600.perlbench_s-210B.champsimtrace.xz` and checked the branch predictors accuracy, MPKI and IPC for 64KB hardware storage. So, to check for all five predictors ie, Bimodal, Gshare, Perceptron, Tage and Hybrid, we do the following:

- Update the  `branch_predictor` field in `champsim_config.json` file in ChampSim folder as the desired name of predictor to be evaluated. For eg: to evaluate gshare, set it to "gshare"

- Then run the command:

```
./config.sh champsim_config.json
make
```
This will update the champsim file in the bin bolder to evaluate gshare (or specified predictor).

- Now, check the simulation for each of your traces using the following command:

```
$ bin/champsim --warmup_instructions 200000000 --simulation_instructions 500000000 ~/path/to/traces/cloud9_phase0_core1.trace.xz -c
```
