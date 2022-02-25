import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math 

pdb_sample = pd.read_fwf(
    "PDB_input_5.txt",
    header=None,
    infer_nrows=200,
)
col_names = [
    "header",
    "id",
    "Atom_id",
    "res_id",
    "strand",
    "res_number",
    "x",
    "y",
    "z",
    "q_num1",
    "q_num2",
    "atom_identifier",
]

Num_atom = int(pd.read_fwf('equill.lammpstrj', header=None, delimiter=' ').iloc[3][0])
trajectory_df = pd.read_fwf('equill.lammpstrj', header=None, delimiter=' ')

Snap_shot_number = int(int(trajectory_df.iloc[::-1].iloc[Num_atom + 7][0]) / (int(trajectory_df.iloc[::-1].iloc[Num_atom + 7][0]) - int(trajectory_df.iloc[::-1].iloc[((Num_atom + 7)*2)+2][0])))
Snap_shot_number
overall = pd.DataFrame()
step_counter = 1
for snapshot in range(Snap_shot_number): 
    n = (Num_atom * (step_counter-1))+((step_counter)*9)
    trajectory_data = pd.read_csv("equill.lammpstrj",skiprows=n,header=None,nrows=Num_atom,delimiter=" ", usecols=[3, 4, 5])
    raw_pdb = pd.concat(
    [pdb_sample[[0, 1, 2, 3, 4, 5]], trajectory_data[[3, 4, 5]], pdb_sample[[9, 10, 11]]], axis=1
    )
    raw_pdb.columns = col_names
    CA_list = raw_pdb.query("Atom_id == 'CA'")
    CA = list(zip(CA_list["res_number"].values, CA_list["x"].values,CA_list["y"].values,CA_list["z"].values))
    Dist_list = []
    Data = pd.DataFrame(columns=['res_num_i', 'res_num_j', 'i', 'j', 'dist'])
    for res_num_i, x_i, y_i, z_i  in CA:
        for res_num_j, x_j, y_j, z_j in CA:
            calc_dist = math.dist([x_i, y_i, z_i], [x_j, y_j, z_j])
            if calc_dist < 6:
                Data = Data.append({'res_num_i': res_num_i, 'res_num_j': res_num_j ,
                'i': [x_i, y_i, z_i], 'j': [x_j, y_j, z_j], 'dist': calc_dist},
                ignore_index=True) # ignore index : index haye ghabl ro dar nazar nagir

                overall = overall.append(
                    {'model' : step_counter, 'i': res_num_i, 'j': res_num_j,'distance': calc_dist},
                     ignore_index=True)
                overall.to_csv('plotting_material.csv')
                
    x = Data["res_num_i"]
    y = Data["res_num_j"]
    plt.figure(figsize=(12,8))
    plt.scatter(x,y, marker="^", cmap="viridis", color="k")
    plt.xlabel("Residue index Number", fontdict={'fontname': 'Times'}, size = 15)
    plt.ylabel("Residue index Number", fontdict={'fontname': 'Times'}, size = 15)
    plt.gca().invert_yaxis()
    plt.title("C"+"\u03B1"+"-"+"C"+"\u03B1"+" Protein Contact Map", fontdict={'fontname': 'Times'}, size = 25)
    name =  "Protein Contact Map - " + str(step_counter)
    plt.savefig(name)
    step_counter += 1