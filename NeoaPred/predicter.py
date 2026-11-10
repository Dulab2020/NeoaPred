import NeoaPred.masif_tools.default_config.global_vars
import os
import pandas as pd
from NeoaPred.PepConf.utils.seq2pdb import seq2pdb
from NeoaPred.PepConf.structure import structure_predicter
from NeoaPred.Surface.surface import compute_surface
from NeoaPred.PepFore.foreignness import foreignness_predicter

def neoantigen_predicter(neo_input_file, output_dir, step, trained_model_1, trained_model_2):
    cons_init_pep = False
    pred_struc = False
    comp_mol_surf = False
    pred_fore = False
    if("0" in step):
        cons_init_pep = True
    if("1" in step):
        pred_struc = True
    if("2" in step):
        comp_mol_surf= True
    if("3" in step):
        pred_fore = True
    
    deal_single_pep = False

    df = pd.read_csv(neo_input_file, header=0)
    
    if('Pep' in df.columns):
        if pred_fore:
            print('The columns name of input file must be \'ID,Allele,WT,Mut,...\' or \'ID,Allele,Pep,...\'')
            print('To predict foreignness score, the columns name of input file must be \'ID,Allele,WT,Mut,...\'')
            os._exit(0)
        deal_single_pep = True
    
    init_pep_dir = output_dir + "/" + "InitPep"
    conf_dir = output_dir + "/" + "Structure"
    surf_dir = output_dir + "/" + "Surface"
    for_dir = output_dir + "/" + "Foreignness"

    initpep_list = []
    data_list = []

    for idx, row in df.iterrows():
        ID = str(row["ID"])
        if deal_single_pep:
            ID_pep = ID + "_pep"
        else:
            ID_wt = ID + "_wt"
            ID_mut = ID + "_mut"
        allele = row["Allele"]
        #init_pep
        if deal_single_pep:
            seq_pep = row["Pep"]
            init_pep_pep = init_pep_dir + "/" + ID_pep + ".pdb"
            initpep_list.append([init_pep_pep, seq_pep])
        else:
            seq_wt = row["WT"]
            init_pep_wt = init_pep_dir + "/" + ID_wt + ".pdb"
            initpep_list.append([init_pep_wt, seq_wt])
            seq_mut = row["Mut"]
            init_pep_mut = init_pep_dir + "/" + ID_mut + ".pdb"
            initpep_list.append([init_pep_mut, seq_mut])
        
        #conformation
        if deal_single_pep:
            conf_pep = conf_dir + "/" + ID_pep +"_relaxed.pdb"
            data_list.append([ID_pep, allele, seq_pep, ID])
        else:
            conf_wt = conf_dir + "/" + ID_wt +"_relaxed.pdb"
            data_list.append([ID_wt, allele, seq_wt, ID])
            conf_mut = conf_dir + "/" + ID_mut +"_relaxed.pdb"
            data_list.append([ID_mut, allele, seq_mut, ID])
    
    if(cons_init_pep == True):
        #Construct initialized peptide
        if not os.path.exists(init_pep_dir):
            os.makedirs(init_pep_dir)
        for out_pdb, seq in initpep_list:
            seq2pdb(seq, out_pdb)

    if(pred_struc == True):
        #Generate MHCI peptide complex structure.
        if not os.path.exists(conf_dir):
            os.makedirs(conf_dir)
        input_file = conf_dir + "/" + "sample_info.csv"
        df_info = pd.DataFrame(data_list, columns = ["ID", "Allele", "Peptide_Seq", "SampleID"])
        df_info.to_csv(input_file, index=False)
        output_file = conf_dir + "/" + "MhcPepStruc_pLDDTs.csv"
        
        structure_predicter(
                            input_dir = init_pep_dir,
                            input_file = input_file,
                            output_dir = conf_dir,
                            output_file = output_file,
                            pre_train_model = trained_model_1,
                            MHC_inDir = None
                           )
    
    if(comp_mol_surf == True):
        #Compute molecular surfaces features.
        if not os.path.exists(surf_dir):
            os.makedirs(surf_dir)
        input_file = surf_dir + "/" + "sample_info.csv"
        df_info = pd.DataFrame(data_list, columns = ["ID", "Allele", "Peptide_Seq", "SampleID"])
        df_info.to_csv(input_file, index=False)
        output_file = surf_dir + "/" + "surface.csv"
 
        compute_surface(
                            input_dir = conf_dir,
                            input_file = input_file,
                            output_dir = surf_dir,
                            output_file = output_file,
                            threads = 3,
                       )

    if(pred_fore == True):
        #Predict foreignness.
        if not os.path.exists(for_dir):
            os.makedirs(for_dir)
        for_cache_dir = surf_dir + "/Cache"
        input_file = neo_input_file
        output_file = for_dir + "/" + "MhcPep_foreignness.csv"

        foreignness_predicter(
                            input_dir = for_cache_dir,
                            input_file = input_file,
                            output_dir = for_dir,
                            output_file = output_file,
                            pre_train_model = trained_model_2,
                          )

