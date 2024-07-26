# import packages
from pathlib import Path
import os
import yaml
import argparse
import pickle
import pandas as pd
import scenicplus
import pycisTopic
import pyranges as pr

# define the defaults here
class DEFAULTS():
    config = {
        "compute_topics": True,
        "enhancer_iden": True,
        "pycisTarget": True,
        "scenicplus": True,
        "scenicplus_downstream": True,
        "load_objects": True,
        "data_path": "./data",
        "input_path": "./data/raw",
        "output_path": "./output",
        "tmp_dir_path": "./gpfs/scratch/your_username",
        "cistarget_db_path": ".data/cisTargetDB",
        "rna_h5ad_name": "rna.h5ad",
        "atac_txt_name": "atac.txt",
        "cell_metadata_name": "metadata.txt",
        "pca_name": "pca.txt",
        "grouping_var": "cell_type",
        "nTopics": None,
        "cpus": 4
    }

# set up configuration of pipeline
def set_up_configs():
    parser = argparse.ArgumentParser(description = "SCENIC+ pipeline for internal use.")
    # input yaml file
    parser.add_argument("config_file",  help = "path to config file with all the options specified", const = None, nargs='?')
    args = parser.parse_args()
    if args.config_file != None:
        config_file = yaml.safe_load(open(args.config_file))
        DEFAULTS.config.update(config_file)

    # # test if inside config tmp dir is the default value
    # if config['tmp_dir_path'] == "./gpfs/scratch/your_username":
    #     # Get the Slurm user ID
    #     slurm_user_id = os.getenv('SLURM_JOB_USER')
    #     print(f'Slurm User ID: {slurm_user_id}')
    #     config['tmp_dir_path'] = f'./gpfs/scratch/{slurm_user_id}'

        
    # # commenting out the below lines while I write the actual pipeline
    # # Check which pipelines (compute_topics, enhancer_iden, pycisTarget, or scenicplus) are specified
    # # to True in the config file
    # if config['compute_topics']:
    #     required_options = ["data_path", "input_path", "output_path", "tmp_dir_path",
    #                     "atac_txt_name", "cell_metadata_name"]
    #     # Check if all required options are present
    #     for option in required_options:
    #         if option not in config or not config[option]:
    #             raise ValueError(f"Required option '{option}' is missing in order to run 'compute topics'")

    # if config['enhancer_iden']:
    #     required_options = ["output_path", "tmp_dir_path"]
    #     # Check if all required options are present
    #     for option in required_options:
    #         if option not in config or not config[option]:
    #             raise ValueError(f"Required option '{option}' is missing in order to run 'enhancer identification'")
        
    #     cistopic_obj_path = Path(config['tmp_dir_path'] / 'cistopic_obj.pkl')
    #     if not cistopic_obj_path.exists():
    #         raise ValueError(f"Required file '{cistopic_obj_path}' from pycisTopic is missing in order to run 'enhancer identification'")


    # if config['pycisTarget']:
    #     required_options = ["output_path", "tmp_dir_path", "cistarget_db_path", 
    #                         ]
    #     # Check if all required options are present
    #     for option in required_options:
    #         if option not in config or not config[option]:
    #             raise ValueError(f"Required option '{option}' is missing in order to run 'enhancer identification'")

    #     required_files = []

    # if config['scenicplus']:
    #     # Code for scenicplus pipeline
    #     pass

    return DEFAULTS

# define the saving function to save keystrokes later
def save_pickle(object, path):
    with open(path) as f:
        pickle.dump(object, f)

def compute_topics(config):
    print("Computing Topics") 
    # get the paths from the config file
    atac_path = Path(f'{config["input_path"]}/{config["atac_txt_name"]}')
    cell_metadata_path = Path(f'{config["input_path"]}/{config["cell_metadata_name"]}')
    tmpDir = Path(config["tmp_dir_path"])
    save_dir = Path(f'{config["tmp_dir_path"]}/scATAC')
    # define some standalone variables for ease
    cpus = int(config['cpus'])

    # read in ATAC count matrix
    count_matrix=pd.read_csv(atac_path, 
                            sep='\t', 
                            engine="pyarrow", 
                            header=0, 
                            index_col=0)
    
    # read in cellular meta data
    cell_data=pd.read_csv(cell_metadata_path, sep='\t')

    # create the cistopic object
    cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix)
    cistopic_obj.add_cell_data(cell_data)

    # topic modelling
    # very computatioanlly intensive!
    # purpose: 
    # 1. To find sets of co-accessible regions (topics), this will be used downstream 
    #   as candidate enhancers (together with Differentially Accessible Regions (DARs)).
    # 2. To impute dropouts. since data is very sparse...
    # Let's generate models with increase ntopics since we are not sure what is correct n (much like PCA)
    # afterwards, choose the model with the "optimal" number of topics (automatically performed in the function)
    models=run_cgs_models(cistopic_obj,
                        n_topics=[2,5,10,15,20,25,30,35,40,45],
                        n_cpu= cpus, # selecting computers here
                        n_iter=500,
                        random_state=555,
                        alpha=50,
                        alpha_by_topic=True,
                        eta=0.1,
                        eta_by_topic=False,
                        save_path=None,
                        _temp_dir = tmpDir)

    # # save results
    # if not save_dir.exists():
    #     save_dir.mkdir()

    # save the cistopic object
    save_pickle(cistopic_obj, save_dir / 'cistopic_obj.pkl')

    print("Finished computing topics, evaluating models")
    # 4 quality control metrics to select the correct model    

    # # uncomment if loading in models from abaove (straight through is commented out)
    # models = pickle.load(open(os.path.join(tmpDir, 'scATAC/models/control_skin_models_500_iter_LDA.pkl'), 'rb'))
    # cistopic_obj = pickle.load(open(os.path.join(tmpDir, 'cistopic_obj.pkl'), 'rb'))
    model = evaluate_models(models,
                           return_model=True,
                           metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                           plot = True,
                           save = save_dir / 'models/evaluated_models.png')
    cistopic_obj.add_LDA_model(model)
    # save the object
    save_pickle(cistopic_obj, save_dir / 'cistopic_obj.pkl')

    print("Finished evaluating models, saving results")
    return(cistopic_obj)

def indentify_enhancers(config, cistopic_obj):
    save_dir = Path(f'{config["tmp_dir_path"]}/scATAC')
    grouping_var = config['grouping_var']

    region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
    region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

    imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
    variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
    markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable= grouping_var, var_features=variable_regions, split_pattern = '-')

    # save the results
    save_pickle(region_bin_topics_otsu, save_dir / 'region_bin_topics_otsu.pkl')
    save_pickle(region_bin_topics_top3k, save_dir / 'region_bin_topics_top3k.pkl')
    save_pickle(markers_dict, save_dir / 'markers_dict.pkl')

def run_pycisTarget(config, region_bin_topics_otsu, region_bin_topics_top3k, markers_dict):
    from scenicplus.wrappers.run_pycistarget import run_pycistarget
    save_dir = Path(f'{config["tmp_dir_path"]}/scATAC')
    cistarget_db_path = Path(config['cistarget_db_path'])
    tmpDir = Path(config["tmp_dir_path"])

    region_sets = {}
    region_sets['topics_otsu'] = {}
    region_sets['topics_top_3'] = {}
    region_sets['DARs'] = {}
    for topic in region_bin_topics_otsu.keys():
        regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    for topic in region_bin_topics_top3k.keys():
        regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    for DAR in markers_dict.keys():
        regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

    rankings_db = Path(cistarget_db_path) /  'motifs-.regions_vs_motifs.rankings.feather'
    scores_db = Path(cistarget_db_path) /  'motifs-.regions_vs_motifs.scores.feather'
    motif_annotation = Path(cistarget_db_path) /  'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

    annot_og = pd.read_csv(os.path.join(db_fpath, "homo-sapien-annot.txt"), delimiter = "\t")
    annot_og.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    run_pycistarget(
        region_sets = region_sets,
        # species = 'homo_sapiens',
        save_path = save_dir,
        ctx_db_path = rankings_db,
        dem_db_path = scores_db,
        path_to_motif_annotations = motif_annotation,
        run_without_promoters = True,
        n_cpu = cpus,
        _temp_dir = os.path.join(tmpDir, 'ray_spill'),
        annotation_version = 'v10nr_clust',
        species = "custom", # species has to equal custom for our own annotations (I think)
        custom_annot=annot_og
        )

def run_scenicplus(config, cistopic_obj,adata,pca_df, menr, tf_file):
    save_dir = Path(f'{config["tmp_dir_path"]}/scATAC')
    cpus = int(config['cpus'])

    scplus_obj = create_SCENICPLUS_object(
        GEX_anndata = adata,
        cisTopic_obj = cistopic_obj,
        menr = menr,
        bc_transform_func = lambda x: f'{x}___cisTopic' #function to convert scATAC-seq barcodes to scRNA-seq ones
    )
    scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())

    #only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
    scplus_obj.dr_cell['GEX_X_pca'] = pca_df
    scplus_obj.dr_cell['GEX_rep'] = pca_df

    # generate cistromes
    merge_cistromes(scplus_obj)

    # Identify the biomart host with the most overlap (done in jupyter notebook but defined here)
    biomart_host = "http://sep2019.archive.ensembl.org/"

    # run enhancer to gene enrichments
    # define the search space around genes: default +- 150kb around TSS
    get_search_space(scplus_obj,
                 biomart_host = biomart_host,
                 species = 'hsapiens',
                 assembly = 'hg38',
                 upstream = [1000, 150000],
                 downstream = [1000, 150000])
    # build enhancer to gene models
    calculate_regions_to_genes_relationships(scplus_obj,
                    ray_n_cpu = cpus,
                    _temp_dir = os.path.join(tmpDir, 'ray_spill'),
                    importance_scoring_method = 'GBM',
                    importance_scoring_kwargs = GBM_KWARGS)

    # infer TF to gene relationships
    calculate_TFs_to_genes_relationships(scplus_obj,
                    tf_file = tf_file,
                    ray_n_cpu = cpus,
                    method = 'GBM',
                    _temp_dir = os.path.join(tmpDir, 'ray_spill'),
                    key= 'TF2G_adj')
    
    # build eGRNs
    build_grn(scplus_obj,
         min_target_genes = 5, 
         adj_pval_thr = 1,
         min_regions_per_gene = 0, 
         quantiles = (0.80, 0.85, 0.90, 0.95),
         top_n_regionTogenes_per_gene = (5, 10, 15, 20),
         top_n_regionTogenes_per_region = (),
         binarize_using_basc = True,
         rho_dichotomize_tf2g = True,
         rho_dichotomize_r2g = True,
         rho_dichotomize_eregulon = True,
         rho_threshold = 0.05,
         keep_extended_motif_annot = True,
         merge_eRegulons = True,
         order_regions_to_genes_by = 'importance',
         order_TFs_to_genes_by = 'importance',
         key_added = 'eRegulons_importance',
         cistromes_key = 'Unfiltered',
         disable_tqdm = False, #If running in notebook, set to True
         ray_n_cpu = cpus,
         _temp_dir = os.path.join(tmpDir, 'ray_spill'))

    # save the results
    save_pickle(scplus_obj, save_dir / 'scenicplus_obj.pkl')
    return(scenicplus_obj)

def run_downstream_scenicplus(config, scplus_obj):
    
    grouping_var = config['grouping_var']
    # format the GRN data frame
    format_egrns(scplus_obj, eregulons_key = 'eRegulons_importance', TF2G_key = 'TF2G_adj', key_added = 'eRegulon_metadata')

    # doing cell type specific analysis
    get_eRegulons_as_signatures(scplus_obj, eRegulon_metadata_key='eRegulon_metadata', key_added='eRegulon_signatures')

    region_ranking = make_rankings(scplus_obj, target='region')
    # Score region regulons
    score_eRegulons(scplus_obj,
                    ranking = region_ranking,
                    eRegulon_signatures_key = 'eRegulon_signatures',
                    key_added = 'eRegulon_AUC',
                    enrichment_type= 'region',
                    auc_threshold = 0.05,
                    normalize = False,
                    n_cpu = 1)
    gene_ranking = make_rankings(scplus_obj, target='gene')
    # Score gene regulons
    score_eRegulons(scplus_obj,
                    gene_ranking,
                    eRegulon_signatures_key = 'eRegulon_signatures',
                    key_added = 'eRegulon_AUC',
                    enrichment_type = 'gene',
                    auc_threshold = 0.05,
                    normalize= False,
                    n_cpu = 1)
    
    # building TF--eGRN relationships
    generate_pseudobulks(scplus_obj,
                         variable = f'GEX_{grouping_var}',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Gene_based',
                         nr_cells = 5,
                         nr_pseudobulks = 100,
                         seed=555)
    generate_pseudobulks(scplus_obj,
                            variable = f'GEX_{grouping_var}',
                            auc_key = 'eRegulon_AUC',
                            signature_key = 'Region_based',
                            nr_cells = 5,
                            nr_pseudobulks = 100,
                            seed=555)
    TF_cistrome_correlation(scplus_obj,
                        variable = f'GEX_{grouping_var}',
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Gene_based',
                        out_key = f'GEX_{grouping_var}_eGRN_gene_based')
    TF_cistrome_correlation(scplus_obj,
                            variable = f'GEX_{grouping_var}',
                            auc_key = 'eRegulon_AUC',
                            signature_key = 'Region_based',
                            out_key = f'GEX_{grouping_var}_eGRN_region_based')

    # Filtering the eGRNs computed
    df1 = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
    df2 = scplus_obj.uns['eRegulon_AUC']['Region_based'].copy()
    df1.columns = [x.split('_(')[0] for x in df1.columns]
    df2.columns = [x.split('_(')[0] for x in df2.columns]
    correlations = df1.corrwith(df2, axis = 0)
    correlations = correlations[abs(correlations) > 0.6]
    # Kepp only R2G +
    keep = [x for x in correlations.index if '+_+' in x] + [x for x in correlations.index if '-_+' in x]
    # Keep extended if not direct
    extended = [x for x in keep if 'extended' in x]
    direct = [x for x in keep if not 'extended' in x]
    keep_extended = [x for x in extended if not x.replace('extended_', '') in direct]
    keep = direct + keep_extended
    # Keep regulons with more than 10 genes
    keep_gene = [x for x in scplus_obj.uns['eRegulon_AUC']['Gene_based'].columns if x.split('_(')[0] in keep]
    keep_gene = [x for x in keep_gene if (int(x.split('_(')[1].replace('g)', '')) > 10)]
    keep_all = [x.split('_(')[0] for x in keep_gene]
    keep_region = [x for x in scplus_obj.uns['eRegulon_AUC']['Region_based'].columns if x.split('_(')[0] in keep]
    scplus_obj.uns['selected_eRegulons'] = {}
    scplus_obj.uns['selected_eRegulons']['Gene_based'] = keep_gene
    scplus_obj.uns['selected_eRegulons']['Region_based'] = keep_region

    # binarize eGRNS
    binarize_AUC(scplus_obj,
             auc_key='eRegulon_AUC',
             out_key='eRegulon_AUC_thresholds',
             signature_keys=['Gene_based', 'Region_based'],
             n_cpu=cpus)

    # Run Dimensionlaioty Reduction on the eGRNs
    run_eRegulons_umap(scplus_obj,
                   scale=True, signature_keys=['Gene_based', 'Region_based'], selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'])
    run_eRegulons_tsne(scplus_obj,
                   scale=True, signature_keys=['Gene_based', 'Region_based'], selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'])
    
    # specificty scores
    regulon_specificity_scores(scplus_obj,
                         f'GEX_{grouping_var}',
                         signature_keys=['Gene_based'],
                         selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'],
                         out_key_suffix='_gene_based',
                         scale=False)
    
    # save the results
    save_pickle(scplus_obj, save_dir / 'scenicplus_obj_downstream.pkl')

if __name__ == '__main__':

    # import configs
    DEFAULTS=set_up_configs()
    print("Configuration:")
    for key, value in DEFAULTS.config.items():
        print(f'     {key}: {value}')

    # optionally run the pipelines specified in the config file
    if DEFAULTS.config['compute_topics']:
        cistopic_obj=compute_topics(DEFAULTS.config)

    if DEFAULTS.config['enhancer_iden']:
        if DEFAULTS.config['load_objects']:
            save_dir = Path(f'{DEFAULTS.config["tmp_dir_path"]}/scATAC')
            cistopic_obj = pickle.load(open(save_dir / 'cistopic_obj.pkl', 'rb'))
        
        indentify_enhancers(DEFAULTS.config, cistopic_obj)

    if DEFAULTS.config['pycisTarget']:
        save_dir = Path(f'{config["tmp_dir_path"]}/scATAC')
        if DEFAULTS.config['load_objects']:
            # load the objects
            markers_dict = pickle.load(open(save_dir / 'markers_dict.pkl', 'rb'))
            region_bin_topics_otsu = pickle.load(open(save_dir / 'region_bin_topics_otsu.pkl', 'rb'))
            region_bin_topics_top3k = pickle.load(open(save_dir / 'region_bin_topics_top3k.pkl', 'rb'))

        pycisTarget_obj=run_pycisTarget(config, region_bin_topics_otsu, region_bin_topics_top3k, markers_dict)
    
    if DEFAULTS.config['scenicplus']:
        if DEFAULTS.config['load_objects']:
            save_dir = Path(f'{config["tmp_dir_path"]}/scATAC')
            cistopic_obj = pickle.load(open(save_dir / 'cistopic_obj.pkl', 'rb'))
            adata = sc.read_h5ad(os.path.join(f'data/scenicplus/%DEFAULTS.config["db_prefix"].h5ad'))
            menr = pickle.load(open(save_dir / 'menr.pkl', 'rb'))
            pca_df = pd.read_csv('data/scenicplus/pca.txt', sep='\t', engine="pyarrow", header=0)
            tf_file='data/scenicplus/utoronto_human_tfs_v_1.01.txt'
            
        scenicplus_obj= run_scenicplus(DEFAULTS.config,
                        cistopic_obj,
                        adata,
                        pca_df,
                        menr, 
                        tf_file)
    
    if DEFAULTS.config['scenicplus_downstream']:
        if DEFAULTS.config['load_objects']:
            save_dir = Path(f'{config["tmp_dir_path"]}/scATAC')
            scplus_obj = pickle.load(open(save_dir / 'scenicplus_obj.pkl', 'rb'))
        run_downstream_scenicplus(DEFAULTS.config, scenicplus_obj)
    