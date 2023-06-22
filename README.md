# modelling2_gif folder

## Non Columnar version
- No cols response to other stim than AB - `no_cols_res_to_A.m` `no_cols_res_to_B.m` `no_cols_res_to_BA.m`
- Variability in all weights -  `a1_no_cols_full_var.m` , `a1_no_cols_full_var_long.m`
- Variability in E-E weights - `a1_no_cols_ee_w_var.m`, `a1_no_cols_ee_w_var_long.m`
- Normal version: `a1_no_cols.m` , `a1_no_cols_long.m`
- Variability in all types of weights - `a1_no_cols_full_var.m` `a1_no_cols_full_var_long.m`
- No cols - Thalamic cols inc to 25 cols - `a1_no_cols_thalamic_many_cols.m`, `a1_no_cols_thalamic_many_cols_long.m`


# As of 22 June, 2023
Files to run the simulation for AB
- a1_no_cols_thalamic_many_cols.m
- a1_no_cols_thalamic_many_cols_long.m

Files to run AB trained network for A,B,AB
- no_cols_res_to_A.m
- no_cols_res_to_B.m
- no_cols_res_to_BA.m
Note: Change `weight_matrix = load(starting_batch_path).weight_matrix;` for initial and final response. 

Files to run for final analysis btn A,B,AB,BA:-
- no_col_analysis.m