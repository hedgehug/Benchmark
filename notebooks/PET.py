import statsmodels.stats.multitest as multitest
import scipy.stats as stat
import sys, random, math
import numpy as np
import os
import pandas as pd
import gseapy as gp


random.seed(42)


def extract_deg_rnk(file_name, deg_num=200, ascending=True):
    df = pd.read_csv(file_name, sep='\t', header=None)
    df = df.sort_values(by=1, ascending=ascending)
    top_genes = df[0][:deg_num].to_numpy()
    return top_genes


def read_gmt(file_name):
    file = open(file_name, 'r')
    info = file.readlines()
    pathway_dict = dict()
    for line in info:
        line = line.strip().split('\t')
        pathway_dict[line[0]] = line[2:]
    file.close()
    return pathway_dict


def run_ora(pathway_dict, deg_dict, gene_universe_num, out_dir):
    for item in deg_dict.keys():
        print('Running ORA for genes in ', item)
        run_fisher_test(pathway_dict, deg_dict[item], gene_universe_num, out_dir, key=item)


def run_fisher_test(pathway_dict, gene_set, gene_universe_num, out_dir, key):
    gene_set = set(gene_set)
    pval_list = []
    intersection_list = []
    pathway_list = list(pathway_dict.keys())
    # print('Start performing fisher test!')
    for pathway in pathway_list:
        # print(pathway)
        intersection_num, pval = run_fisher(gene_set, pathway_dict[pathway], gene_universe_num)
        pval_list.append(pval)
        intersection_list.append(intersection_num)
    fdr_list = multitest.fdrcorrection(pval_list, is_sorted=False)[-1]
    # sort by p-value
    pval_list, fdr_list, pathway_list, intersection_list = \
        (list(t) for t in zip(*sorted(zip(pval_list, fdr_list, pathway_list, intersection_list))))

    out_file = open(out_dir+'/'+key+'.ora_result.txt', 'w')
    out_file.write('Pathway\t#Overlap Gene\tp-value\tFDR\tRank\n')
    rank = 1
    for idx in range(len(pathway_list)):
        out_file.write('\t'.join([pathway_list[idx], str(intersection_list[idx]),
                                  str(pval_list[idx]), str(fdr_list[idx]), str(rank)])+'\n')
        rank += 1
    out_file.close()
    print('Results written to ', out_dir+'/'+key+'.ora_result.txt')
    print('*'*20)


def run_enrichr(pathway_dict, deg_dict, gene_universe, permutation_num, permutation_file_name, out_dir):
    for item in deg_dict.keys():
        print('Running Enrichr for genes in ', item)
        run_enrichr_test(pathway_dict=pathway_dict, gene_set=deg_dict[item],
            gene_universe=gene_universe, permutation_num=permutation_num,
            permutation_file_name=permutation_file_name, out_dir=out_dir, key=item)


def run_enrichr_test(pathway_dict, gene_set, gene_universe, permutation_num,
        permutation_file_name, out_dir, key):
    if not os.path.exists(permutation_file_name):
        print(permutation_file_name, 'does not exist, initializing permutation')
        generate_permutation(gene_universe, len(gene_set),
                             permutation_num, pathway_dict, permutation_file_name)

        print('Permutation result written to ', permutation_file_name)
        pathway_mean_rank_dict, pathway_std_rank_dict = read_permutation_file(permutation_file_name)
    else:
        print(permutation_file_name, 'already exists, reading file now...')
        pathway_mean_rank_dict, pathway_std_rank_dict = read_permutation_file(permutation_file_name)

    # run fisher exact test
    all_fisher_result = []
    pval_dict = dict()
    pathway_names = list(pathway_dict.keys())
    for pathway in pathway_names:
        intersection, fisher_pval = run_fisher(set(gene_set) , pathway_dict[pathway], len(gene_universe))
        pval_dict[pathway] = fisher_pval
        all_fisher_result.append([pathway, intersection, fisher_pval])

    # sort by p-value
    all_fisher_result.sort(key=lambda x: x[2])

    # compute the original rank
    rank = 1
    rank_dict = dict()
    for item in all_fisher_result:
        rank_dict[item[0]] = rank
        rank += 1

    # compute z-score
    combined_score_dict = dict()
    combined_score_list = []
    for pathway in pathway_names:
        z_score = (rank_dict[pathway] - pathway_mean_rank_dict[pathway]) / pathway_std_rank_dict[pathway]
        combined_score_dict[pathway] = abs(z_score * np.log10(pval_dict[pathway]))
        combined_score_list.append(combined_score_dict[pathway])

    # sort combined score, descending
    combined_score_list, pathway_names = (list(t) for t in zip(*sorted(
        zip(combined_score_list, pathway_names), reverse=True)))

    # fdr_list = multitest.fdrcorrection(pval_list, is_sorted=False)[-1]

    out_file = open(out_dir+'/'+key+'.enrichr_result.txt', 'w')
    out_file.write('Pathway\tp-value\tCombined score\tRank\n')
    rank = 1
    for idx in range(len(pathway_names)):
        out_file.write('\t'.join([pathway_names[idx], str(pval_dict[pathway_names[idx]]),
                                  str(combined_score_list[idx]), str(rank)]) + '\n')
        rank +=1
    out_file.close()
    print('Results written to ', out_dir+'/'+key+'.enrichr_result.txt')
    print('*' * 20)


def generate_permutation(gene_universe, gene_num, permutation_num, pathway_dict, permutation_file_name):
    all_random_set = []
    # sample N genes for perm_num rounds
    for idx in range(0, permutation_num):
        tmp_geneset = random.sample(gene_universe, gene_num)
        all_random_set.append(set(tmp_geneset))

    # pathway: [all rank in all permutations]
    permutation_rank_dict = dict()
    for idx in range(0, permutation_num):
        # print('Permutation ' + str(idx))
        if idx%50 == 0:
            print('Permutation ' + str(idx))
        tmp_rank_list = []
        for pathway in list(pathway_dict.keys()):
            intersection, fisher_pval = run_fisher(all_random_set[idx], pathway_dict[pathway], len(gene_universe))
            tmp_rank_list.append([pathway, intersection, fisher_pval])
        # sort by p-value
        tmp_rank_list.sort(key=lambda x: x[2])
        rank = 1
        for item in tmp_rank_list:
            try:
                permutation_rank_dict[item[0]].append(rank)
            except:
                permutation_rank_dict[item[0]] = []
                permutation_rank_dict[item[0]].append(rank)
            rank += 1

    out_file = open(permutation_file_name, 'w')
    for pathway in permutation_rank_dict.keys():
        tmp_list = permutation_rank_dict[pathway]
        tmp_list.sort()
        out_file.write(pathway + '\t' + '\t'.join([str(i) for i in tmp_list]) + '\n')
    out_file.close()


def run_fisher(geneset_1, geneset_2, gene_universe_num):
    # input gene list
    geneset_num_1 = len(geneset_1)
    # pathway
    geneset_num_2 = len(geneset_2)
    intersection = len(geneset_1.intersection(geneset_2))
    _, fisher_pval = stat.fisher_exact([[intersection, geneset_num_1 - intersection], [geneset_num_2 - intersection,
                                                                                        gene_universe_num - geneset_num_2 - geneset_num_1 + intersection]],
                                        alternative='greater')
    return intersection, fisher_pval


def read_permutation_file(permutation_file):
    file = open(permutation_file, 'r')
    file_info = file.readlines()
    pathway_mean_rank_dict = dict()
    pathway_std_rank_dict = dict()
    for line in file_info:
        line = line.strip().split('\t')
        pathway_mean_rank_dict[line[0]] = np.mean(np.array(line[1:]).astype(int))
        pathway_std_rank_dict[line[0]] = np.std(np.array(line[1:]).astype(int))
    file.close()
    return pathway_mean_rank_dict, pathway_std_rank_dict



def extract_ora_rank(ora_result_file):
    print('Calculating ora rank...')
    file = open(ora_result_file, 'r')
    file_info = file.readlines()
    pathway_pval_dict = dict()
    pathway_rank_dict = dict()
    pathway_list = []
    pval_list = []
    for line in file_info[1:]:
        # Pathway	#Overlap Gene	p-value	FDR	Rank
        line = line.strip().split('\t')
        # pathway, intersection number/combined score, p-value
        pathway_pval_dict[line[0].upper()] = float(line[2])
        pathway_list.append(line[0].upper())
        pval_list.append(float(line[2]))
    # sort by p-value, ascending
    # for ora test, there might be lots of tie of pathways
    # here we want to keep the pathways with same p-value/combined score the same rank, for later average calculation
    pval_ranks = stat.rankdata(pval_list, method='average')
    for idx in range(len(pval_ranks)):
        pathway_rank_dict[pathway_list[idx]] = pval_ranks[idx]
    return pathway_rank_dict, pathway_pval_dict


def extract_enrichr_rank(enrichr_result_file):
    print('Calculating enrichr rank...')
    file = open(enrichr_result_file, 'r')
    file_info = file.readlines()
    pathway_pval_dict = dict()
    pathway_rank_dict = dict()
    pathway_list = []
    combined_score_list = []
    for line in file_info[1:]:
        # Pathway	p-value	Combined score	Rank
        line = line.strip().split('\t')
        pathway_pval_dict[line[0].upper()] = float(line[1])
        pathway_list.append(line[0].upper())
        combined_score_list.append(float(line[2]))

    # sort by combined, descending, and keep the ties
    combined_score_ranks = stat.rankdata(combined_score_list, method='average')
    combined_score_ranks = len(combined_score_ranks)-combined_score_ranks+1
    for idx in range(len(combined_score_ranks)):
        pathway_rank_dict[pathway_list[idx]] = combined_score_ranks[idx]
    return pathway_rank_dict, pathway_pval_dict


def extract_gsea_rank(gsea_result_file, gsea_label='up'):
    print('Calculating GSEA rank...')
    pathway_pval_dict = dict()
    pathway_rank_dict = dict()
    gsea_df = pd.read_csv(gsea_result_file, sep=',', header=0, index_col=0)
    if gsea_label == 'up':
        # positive NES first
        gsea_df = gsea_df.sort_values(by='NES', ascending=False)
    else:
        # negative NES first
        gsea_df = gsea_df.sort_values(by='NES', ascending=True)

    for idx in range(gsea_df.shape[0]):
        pathway_rank_dict[gsea_df.iloc[idx][0]] = idx+1
        pathway_pval_dict[gsea_df.iloc[idx][0]] = gsea_df.iloc[idx]['NOM p-val']
    return pathway_rank_dict, pathway_pval_dict


def run(ora_result_file, enrichr_result_file, gsea_result_file, out_dir, file_prefix, pathway_dict, direction):
    # pathway: value
    rank_dict = dict()
    pval_dict = dict()

    rank_dict['ora'], pval_dict['ora'] = extract_ora_rank(ora_result_file)
    rank_dict['Enrichr'], pval_dict['Enrichr'] = extract_enrichr_rank(enrichr_result_file)
    rank_dict['GSEA'], pval_dict['GSEA'] = extract_gsea_rank(gsea_result_file, gsea_label=direction)

    methods = list(rank_dict.keys())
    methods.sort()

    print('Ensembling results...')
    pathway_list = list(pathway_dict.keys())
    results = []
    for p in pathway_list:
        tmp_rank_list = []
        tmp_pval_list = []
        for method in methods:
            try:
                tmp_rank_list.append(rank_dict[method][p])
                tmp_pval_list.append(pval_dict[method][p])
            except:
                print(p, 'is missing in', method, '. Skipped.')
                break
        if len(tmp_rank_list) == len(tmp_pval_list) == len(methods):
            # compute average rank
            avg_rank = np.mean(tmp_rank_list)
            # compute combined p-value
            combined_pval = stat.combine_pvalues(tmp_pval_list, method='stouffer')[-1]
            if np.isnan(combined_pval):
                # certain p-value lists cannot be combined, e.g. [1.0, 1.0, 0.0], assign as 1.0 for FDR correction
                combined_pval = 1.0
            results.append([p, avg_rank, combined_pval])
        else:
            continue

    # sort based on average rank
    results.sort(key=lambda x: x[1])

    # calculate FDR
    pval_list = np.array(results)[:, -1].astype(float)
    fdr_list = multitest.fdrcorrection(pval_list, is_sorted=False)[-1]

    # write to results
    out_file = open(out_dir+'/'+file_prefix+'.PET.txt', 'w')
    # write header
    out_file.write('Pathway\tPET rank\tPET p-value\tPET FDR\tAverage rank')
    for method in methods:
        out_file.write('\t'+method+' rank\t'+method+' p-value')
    out_file.write('\n')
    for idx in range(len(results)):
        out_file.write('\t'.join([results[idx][0], str(idx+1),
                                  str(results[idx][2]), str(fdr_list[idx]), str(results[idx][1])]))
        for method in methods:
            out_file.write('\t' +str(rank_dict[method][results[idx][0]])+'\t'+str(pval_dict[method][results[idx][0]]))
        out_file.write('\n')
    out_file.close()

    print('Results written to ', out_dir+'/'+file_prefix+'.PET.txt')


def combine_gsea_cmd_results(gsea_dir, gsea_label, key):
    # GSEA results are stored in two tsv/csv files
    gsea_report_prefix = 'gsea_report_for_'
    # gsea_report_suffix = 'tsv'

    gsea_report_list = []
    all_file = os.listdir(gsea_dir)
    all_file.sort()
    for f in all_file:
        if f.startswith(gsea_report_prefix) and not f.endswith('html'):
            gsea_report_list.append(gsea_dir + f)

    if len(gsea_report_list) != 2:
        print('Did not find the expected 2 GSEA reports. Please check GSEA run. Stopping now.')
        return 0, 0

    # for the report of expected direction, simply take the rank, for the other report, do the reverse way
    if gsea_report_list[0].__contains__(gsea_label):
        pos_report = gsea_report_list[0]
        neg_report = gsea_report_list[1]
    else:
        pos_report = gsea_report_list[1]
        neg_report = gsea_report_list[0]

    df_1 = pd.read_csv(pos_report, sep='\t')
    df_2 = pd.read_csv(neg_report, sep='\t')
    # for the neg one, correct p-value to 1
    df_2['NOM p-val'] = [1] * df_2.shape[0]
    new_df = pd.concat([df_1, df_2])
    new_df.insert(0, "Name", ['Gsea']*new_df.shape[0])
    # write to a tmp file, delete later
    new_df.to_csv(gsea_dir+'/tmp_result.'+key+'.txt', sep=',', index=False)
    return gsea_dir+'/tmp_result.'+key+'.txt'


def run_PET(ora_result_file, enrichr_result_file, gsea_result_file, pathway_dict,
            out_dir, gsea_pos_label='pos', file_prefix=None):
    # identify files relevant
    files = dict()
    files['up'] = dict()
    files['up']['ora'] = ora_result_file
    files['up']['enrichr'] = enrichr_result_file

    tmp_file = []
    gsea_path = gsea_result_file
    if os.path.isdir(gsea_path):
        print('GSEA path was provided as directory, command line version used.')
        # combine GSEA results, make it GSEAPY output like
        if gsea_pos_label == 'pos':
            files['up']['gsea'] = combine_gsea_cmd_results(gsea_dir=gsea_path, gsea_label=gsea_pos_label, key='up')
            files['down']['gsea'] = combine_gsea_cmd_results(gsea_dir=gsea_path, gsea_label='neg', key='down')
        elif gsea_pos_label == 'neg':
            files['up']['gsea'] = combine_gsea_cmd_results(gsea_dir=gsea_path, gsea_label=gsea_pos_label, key='up')
            files['down']['gsea'] = combine_gsea_cmd_results(gsea_dir=gsea_path, gsea_label='pos', key='down')
        else:
            print('Please provide a valid direction key to identify GSEA report. Options: pos, neg. '
                  'The key should indicate up-regulation in the condition.')
        tmp_file.append(files['up']['gsea'])
        tmp_file.append(files['down']['gsea'])
    else:
        print('GSEA path was provided as file, GSEAPY used.')
        # correct the p-value and delete the tmp files later
        tmp_df = pd.read_csv(gsea_path, sep=',', header=0, index_col=0)
        tmp_df.loc[tmp_df['NES'] <= 0, 'NOM p-val'] = 1
        tmp_df.to_csv(gsea_path.replace('.csv', '.up.tmp.csv'), sep=',', index=True)
        # GSEA might have NA p-value, replace with 1
        tmp_df['NOM p-val'].fillna(1, inplace=True)
        files['up']['gsea'] = gsea_path.replace('.csv', '.up.tmp.csv')
        tmp_file.append(files['up']['gsea'])
        tmp_df = pd.read_csv(gsea_path, sep=',', header=0, index_col=0)
        tmp_df.loc[tmp_df['NES'] >= 0, 'NOM p-val'] = 1
        # GSEA might have NA p-value, replace with 1
        tmp_df['NOM p-val'].fillna(1, inplace=True)
        # tmp_df.to_csv(gsea_path.replace('.csv', '.down.tmp.csv'), sep=',', index=True)
        # files['down']['gsea'] = gsea_path.replace('.csv', '.down.tmp.csv')
        # tmp_file.append(files['down']['gsea'])

    # combine the results
    for key in files.keys():
        run(ora_result_file=files[key]['ora'], enrichr_result_file= files[key]['enrichr'],
            gsea_result_file=files[key]['gsea'], out_dir=out_dir,
            pathway_dict=pathway_dict, direction=key, file_prefix=file_prefix)

    # remove the tmp files created for GSEA
    for f in tmp_file:
        os.remove(f)


def run_GSEA(prerank_file_path, pathway_file, out_dir, gsea_cli_path=None, method='Python',
             thread_num=5, min_size=15, max_size=500, permutation_num=1000, seed=42,
             no_plot=True, gsea_out_label='GSEA_result'):
    if method == 'Python':
        run_GSEA_python(prerank_file_path = prerank_file_path,
                        pathway_file = pathway_file, out_dir=out_dir,
                        no_plot= no_plot, min_size=min_size, max_size=max_size,
                        thread_num = thread_num, permutation_num=permutation_num,
                        seed=seed)
    elif method == 'cli':

        if gsea_cli_path == None:
            raise Exception('Please specify file path of GSEA gsea-cli.sh. Otherwise, please keep method=\'Python\'.')
        if not os.path.exists(gsea_cli_path):
            raise Exception(gsea_cli_path+' file does not exist!')

        # call GSEA cmd script
        run_GSEA_cmd(cli_file_path=gsea_cli_path, prerank_file_path = prerank_file_path,
                     pathway_file = pathway_file, out_dir=out_dir,
                     plot= no_plot, min_size=min_size, max_size=max_size,
                     permutation_num=permutation_num, seed=seed, gsea_out_label=gsea_out_label)
    return


def run_GSEA_cmd(cli_file_path, prerank_file_path, pathway_file, out_dir, gsea_out_label,
                 plot=False, min_size=15, max_size=500, permutation_num=1000, seed=42):
    """
    ~/GSEA_test/GSEA_cmd/gsea-cli.sh GSEAPreranked -gmx example_new/c2.cp.kegg.v2023.1.Hs.symbols.gmt
    -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk example_new/cond1.vs.cond2.rnk
    -scoring_scheme weighted -rpt_label example_test   -create_svgs false -include_only_symbols true
     -make_sets true -plot_top_x 5 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false
     -out example_new/
    :return:
    """
    print('Start running GSEA!')
    print('Commandas to run GSEA: '+cli_file_path+' GSEAPreranked -gmx '+pathway_file+
              ' -collapse No_Collapse -mode Max_probe -norm meandiv -nperm '+str(permutation_num)
              +' -rnk '+prerank_file_path+' -scoring_scheme weighted -rpt_label '+gsea_out_label
              +' -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 5 -rnd_seed '+str(seed)+
              ' -set_max '+str(max_size)+' -set_min '+str(min_size)+' -zip_report false -out '+out_dir)
    os.system(cli_file_path+' GSEAPreranked -gmx '+pathway_file+
              ' -collapse No_Collapse -mode Max_probe -norm meandiv -nperm '+str(permutation_num)
              +' -rnk '+prerank_file_path+' -scoring_scheme weighted -rpt_label '+gsea_out_label
              +' -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 5 -rnd_seed '+str(seed)+
              ' -set_max '+str(max_size)+' -set_min '+str(min_size)+' -zip_report false -out '+out_dir)


def run_GSEA_python(prerank_file_path, pathway_file, out_dir, no_plot=True, thread_num=5, min_size=15,
                    max_size=500, permutation_num=1000, seed=42):
    # read rank file
    rnk = pd.read_csv(prerank_file_path, header=None, index_col=0, sep="\t")
    # # run prerank
    # # enrichr libraries are supported by prerank module. Just provide the name
    # # use 4 process to acceralate the permutation speed
    pre_res = gp.prerank(rnk=rnk,  # or rnk = rnk,
                         gene_sets=pathway_file,
                         threads=thread_num,
                         min_size=min_size,
                         max_size=max_size,
                         permutation_num=permutation_num,  # reduce number to speed up testing
                         outdir=out_dir,  # don't write to disk
                         seed=seed,
                         no_plot = no_plot,
                         verbose=True)
