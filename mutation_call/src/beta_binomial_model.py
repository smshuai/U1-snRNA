#!/usr/bin/env python
import sys
import numpy as np
import scipy as sp
import pandas as pd
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import betaln
from scipy.stats import combine_pvalues

N_normal = 50 # maximum number of normal samples to use in fitting

def get_args():
    ''' Process argmuents
    '''
    if len(sys.argv) == 5:
        tumour_dp_path = sys.argv[1]
        normal_panel_path = sys.argv[2]
        ref_map_path = sys.argv[3]
        out_path = sys.argv[4]
    elif len(sys.argv)==4:
        tumour_dp_path = sys.argv[1]
        normal_panel_path = sys.argv[2]
        ref_map_path = sys.argv[3]
        out_path = ''
    else:
        sys.stderr.write('Usage: beta_binomial_model.py tumour_depth normal_panel ref_map [out_path]')
        sys.exit(1)
    return tumour_dp_path, normal_panel_path, ref_map_path, out_path


def beta_binomial_pmf(params, k, n):
    ''' PMF of beta-binomial distribution
    '''
    llh = beta_binomial_loglikelihood(params, k, n)
    return np.exp(llh)


def beta_binom_pvalue(params, k, n):
    ''' Calculate Prob(X>=k|params, n, k)
    '''
    tempPV = 0
    for kk in range(int(k), int(np.floor(n + 1))):
        currentValue = beta_binomial_pmf(params, kk, n)
        tempPV = tempPV + currentValue
    return tempPV


def beta_binomial_loglikelihood(params, Ks, Ns):
    """ Calculating log-likelihood of beta-binomial distribution

    Args:
        params (List[float]): the parameter of beta distribution ([alpha, beta])
        Ks (numpy.array([int])): the counts for success
        Ns (numpy.array([int])): the counts of trials
    Note: This function is used to replace the one used in EBfilter (more efficient and compatiable with python3)
    """
    alpha = params[0]
    beta = params[1]
    p1 = np.log(sp.special.comb(Ns, Ks))
    p2 = betaln(Ks + alpha, Ns - Ks + beta)
    p3 = betaln(alpha, beta)
    # log-likelihood
    llh = np.sum(p1 + p2 - p3)
    return llh


def regularized_cost(params, Ks, Ns, reg=0.5):
    ''' Cost function for fitting beta binomial distribution with regularization
    '''
    alpha = params[0]
    beta = params[1]
    llh = beta_binomial_loglikelihood(params, Ks, Ns)
    cost = reg * np.log(alpha + beta) - llh
    return cost


def fit_beta_binomial(Ks, Ns, cost=regularized_cost):
    """ Obtaining maximum likelihood estimator of beta-binomial distribution

    Args:
        Ks (numpy.array([int])): the counts for success
        Ns (numpy.array([int])): the counts of trials
    """
    result = fmin_l_bfgs_b(cost, [20, 20], args = (Ks, Ns), approx_grad = True, bounds = [(0.1, 10000000), (1, 10000000)])
    return result[0]


def calc_bb_pval(Ks, Ns, k, n):
    fit = fit_beta_binomial(Ks, Ns)
    pval = beta_binom_pvalue(fit, k, n)
    return pval


def pval_to_eb(pval):
    if pval < 1e-60:
        EB_score = 60
    elif pval > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = - round(np.log10(pval), 3)
    return EB_score


def test_base(Ks_p, Ns_p, Ks_n, Ns_n, k_pt, n_pt, k_nt, n_nt, k_pn, n_pn, k_nn, n_nn):
    """ Use beta binomial model to call mutation for the base
    Training data based on unpaired normal panels: Ks_p, Ns_p, Ks_n, Ns_n
    Test data from tumour: k_pt, n_pt, k_nt, n_nt
    Test data from paired normal: k_pn, n_pn, k_nn, n_nn
    """
    global N_normal
    # Convert data to int
    arraytoint = lambda x: x.round().astype(np.int)
    Ks_p = arraytoint(Ks_p); Ns_p = arraytoint(Ns_p); Ks_n = arraytoint(Ks_n); Ns_n = arraytoint(Ns_n)
    k_pt = np.int(k_pt); k_nt = np.int(k_nt); n_pt = np.int(n_pt); n_nt = np.int(n_nt)
    k_pn = np.int(k_pn); k_nn = np.int(k_nn); n_pn = np.int(n_pn); n_nn = np.int(n_nn)
    # Filter Ks and Ns by coverage (>=12 reads) and validity (Ks/Ns < 0.5)
    keep_p = np.flatnonzero(np.logical_and(Ns_p>=10, Ks_p / (Ns_p + 0.1) < 0.5))
    keep_n = np.flatnonzero(np.logical_and(Ns_n>=10, Ks_n / (Ns_n + 0.1) < 0.5))
    # Don't need too many training normals
    keep_p = np.random.choice(keep_p, N_normal, replace=False) if keep_p.shape[0] > N_normal else keep_p
    keep_n = np.random.choice(keep_n, N_normal, replace=False) if keep_n.shape[0] > N_normal else keep_n
    assert keep_p.shape[0] >= 30, "<30 useful normal samples for positive strand"
    assert keep_n.shape[0] >= 30, "<30 useful normal samples for negative strand"
    # Avoid that # success > # trials
    if k_pt > n_pt:
        k_pt = n_pt
    if k_nt > n_nt:
        k_nt = n_nt
    if k_pn > n_pn:
        k_pn = n_pn
    if k_nt > n_nt:
        k_nn = n_nn
    # fit model and test
    fitp = fit_beta_binomial(Ks_p[keep_p], Ns_p[keep_p])
    fitn = fit_beta_binomial(Ks_n[keep_n], Ns_n[keep_n])
    pval_pt = beta_binom_pvalue(fitp, k_pt, n_pt)
    pval_nt = beta_binom_pvalue(fitn, k_nt, n_nt)
    pval_pn = beta_binom_pvalue(fitp, k_pn, n_pn)
    pval_nn = beta_binom_pvalue(fitn, k_nn, n_nn)
    cap_pval = lambda p: 1e-60 if p < 1e-60 else p
    pvalt = combine_pvalues([cap_pval(pval_pt), cap_pval(pval_nt)], 'fisher')[1]
    pvaln = combine_pvalues([cap_pval(pval_pn), cap_pval(pval_nn)], 'fisher')[1]
    ebt = pval_to_eb(pvalt); ebn = pval_to_eb(pvaln)
    return ebt, ebn, int(keep_p.shape[0]), int(keep_n.shape[0])


if __name__ == '__main__':
    tumour_dp_path, normal_panel_path, ref_map_path, out_path = get_args()
    tumour_depth = pd.read_table(tumour_dp_path)
    normal_panel = pd.read_table(normal_panel_path)
    ref_map = pd.read_table(ref_map_path)  # reference base at each position
    # Check donor ID
    donor_id = tumour_depth['id'][0]
    if donor_id in normal_panel['id'].values:
        # Get paired normal depth
        normal_depth = normal_panel[normal_panel.id == donor_id]
        normal_panel = normal_panel[normal_panel.id != donor_id]
        # Get output path
        if out_path is '':
            out_path = donor_id + "_U1_GT_res.tsv"
    else:
        sys.stderr.write('Paired normal sample cannot be found for this ID: ' + donor_id)
        sys.exit(1)
    normal_panel = normal_panel.set_index(['chrom', 'pos', 'alt']).sort_index()
    # Add ref to tumour_depth and normal_depth
    tumour_depth = pd.merge(ref_map, tumour_depth, on=('chrom', 'pos'))
    tumour_depth = tumour_depth.loc[tumour_depth.ref != tumour_depth.alt]
    tumour_depth['varDP'] = tumour_depth.varP + tumour_depth.varN
    tumour_depth['totDP'] = tumour_depth.dpP + tumour_depth.dpN
    normal_depth = pd.merge(ref_map, normal_depth, on=('chrom', 'pos'))
    normal_depth = normal_depth.loc[normal_depth.ref != normal_depth.alt]
    normal_depth['varDP'] = normal_depth.varP + normal_depth.varN
    normal_depth['totDP'] = normal_depth.dpP + normal_depth.dpN
    # Merge tumour_depth and normal_depth
    res = pd.merge(tumour_depth, normal_depth, on=('chrom', 'pos', 'ref', 'alt', 'id', 'start', 'end', 'gene', 'strand'), suffixes=('_tumour', '_normal'))
    # Test for tumour and normal sample
    eb_scores_tumour = []
    eb_scores_normal = []
    num_normal_p = []  # number of normal samples used in model for pos strand
    num_normal_n = []  # number of normal samples used in model for neg strand
    for ix, row in res.iterrows():
        # Only test when tumour has variant reads
        if row.varDP_tumour > 0:
            # Get observed number (k_pt: # variant for + strand and tumour)
            k_pt, n_pt, k_nt, n_nt = row[['varP_tumour', 'dpP_tumour', 'varN_tumour', 'dpN_tumour']]
            k_pn, n_pn, k_nn, n_nn = row[['varP_normal', 'dpP_normal', 'varN_normal', 'dpN_normal']]
            # Get normal panel depths
            Ks_p = normal_panel.loc[(row.chrom, row.pos, row.alt), 'varP'].values
            Ns_p = normal_panel.loc[(row.chrom, row.pos, row.alt), 'dpP'].values
            Ks_n = normal_panel.loc[(row.chrom, row.pos, row.alt), 'varN'].values
            Ns_n = normal_panel.loc[(row.chrom, row.pos, row.alt), 'dpN'].values
            # Test tumour and normal with the same model
            ebt, ebn, num_p, num_n = test_base(Ks_p, Ns_p, Ks_n, Ns_n, k_pt, n_pt, k_nt, n_nt, k_pn, n_pn, k_nn, n_nn)
            eb_scores_tumour.append(ebt)
            eb_scores_normal.append(ebn)
            num_normal_p.append(num_p)
            num_normal_n.append(num_n)
        else:
            # No variant read
            eb_scores_tumour.append(0)
            eb_scores_normal.append(0)
            num_normal_p.append(-1)
            num_normal_n.append(-1)
    res['num_p'] = num_normal_p
    res['num_n'] = num_normal_n
    res['EB_tumour'] = eb_scores_tumour
    res['EB_normal'] = eb_scores_normal
    res['EB_delta'] = res.EB_tumour.subtract(res.EB_normal)
    # Assign GT to each pos
    res['GT'] = 'WT'
    # Mut should have large EB in tumour and small EB in normal
    # And the difference in EB between tumour and normal should also be large enough
    Mut = np.logical_and(res.EB_normal<2, res.EB_tumour>4)
    Mut = np.logical_and(Mut, res.EB_delta>4)
    res.loc[Mut, 'GT'] = 'MUT'
    # undet should be
    # 1. not enough depth
    # 2. uncertain number of var reads in normal
    res.loc[res.totDP_tumour < 15, 'GT'] = 'undet'
    res.loc[np.logical_and(Mut, res.varDP_normal>2), 'GT'] = 'undet'
    res.to_csv(out_path, sep='\t', index=False)
