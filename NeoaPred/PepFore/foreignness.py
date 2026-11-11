#Predict Foreignness

import os
import numpy as np
import pandas as pd
import torch
import torch.utils.data as Data
from .dataread import preload, DataSet
from .model import Predict

def evaluate(model, df, device, batch_size, cache_data):
    model = model.to(device)
    model_eval = model.eval()
    torch.manual_seed(99)
    torch.cuda.manual_seed(99)
    eval_ = DataSet(df, cache_data)
    _size = df.shape[0]
    eval_loader = Data.DataLoader(dataset=eval_, batch_size=batch_size, \
                                   shuffle=False, num_workers=2, drop_last=False)

    id_list = []
    eval_list = []
    labels_eval_list = []
    with torch.no_grad():
        for batch_i, (IDs, rho_wt, theta_wt, feat_wt, mask_wt, 
                           rho_mut, theta_mut, feat_mut, mask_mut, 
                           patch_dist, patch_mask, patch_dist_mp, 
                           dist, dist_mask, atom_wt, atom_mut) in enumerate(eval_loader):
            id_list.extend(IDs)
            #wt
            rho_wt = rho_wt.to(device)
            theta_wt = theta_wt.to(device)
            feat_wt = feat_wt.to(device)
            mask_wt = mask_wt.to(device)
            #mut
            rho_mut = rho_mut.to(device)
            theta_mut = theta_mut.to(device)
            feat_mut = feat_mut.to(device)
            mask_mut = mask_mut.to(device)
            #
            patch_mask = patch_mask.to(device)
            patch_dist = patch_dist.to(device)
            patch_dist_mp = patch_dist_mp.to(device)
            #dist
            dist = dist.to(device)
            dist_mask = dist_mask.to(device)
            atom_wt = atom_wt.to(device)
            atom_mut = atom_mut.to(device)
            
            values_out = model(rho_wt, theta_wt, feat_wt, mask_wt, 
                               rho_mut, theta_mut, feat_mut, mask_mut,
                               patch_dist, patch_mask, patch_dist_mp, 
                               dist, dist_mask, atom_wt, atom_mut)

            eval_list.extend(values_out.cpu().numpy())
            labels_eval = values_out.clone()
            labels_eval[torch.where(labels_eval>=0.5)[0]] = 1.0
            labels_eval[torch.where(labels_eval<0.5)[0]] = 0.0
            labels_eval_list.extend(labels_eval.cpu().numpy())
    #df_res = pd.DataFrame([id_list, labels_eval_list, eval_list]).T
    df_res = pd.DataFrame([id_list, eval_list]).T
    df_res.columns = ['SampleID', 'Foreignness_Score']
    return df_res

def foreignness_predicter(
                        input_dir,
                        input_file,
                        output_dir,
                        output_file,
                        pre_train_model,
                        device = "cpu",
                        batch_size = 4,
                        input_df = None,
                      ):
    cache_data = preload(input_dir)
    if input_df is None:
        df = pd.read_csv(input_file, header=0, dtype={"ID": str})
    else:
        df = input_df.copy()
    out_for = output_dir + "/" + "foreignness_result.csv"

    model = Predict(device=device)
    model.load_state_dict(torch.load(pre_train_model, map_location=device))
    df_res = evaluate(model, df, device, batch_size, cache_data)
    df_res.to_csv(out_for, index=False)
    data_pair_list = []
    for sam_ID in df["ID"].unique():
        allele = df[df["ID"]==sam_ID]["Allele"].values[0]
        seq_wt = df[df["ID"]==sam_ID]["WT"].values[0]
        seq_mut = df[df["ID"]==sam_ID]["Mut"].values[0]
        score = df_res[df_res["SampleID"]==sam_ID]["Foreignness_Score"].values[0]
        data_pair_list.append([sam_ID, allele, seq_wt, seq_mut, score])
    df_pair = pd.DataFrame(data_pair_list, columns = ["ID", "Allele", "WT", "Mut", "Foreignness_Score"])
    df_pair.to_csv(output_file, index=False)
