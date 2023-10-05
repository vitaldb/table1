import os
os.environ['R_HOME'] = r'C:\Program Files\R\R-4.3.1'
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stat
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2rpy as np2r

def format_pval(pval):
    '''Returns the formated string of a p-value'''
    if pval < 0.001:
        return '< 0.001'
    return f'{pval:.3f}'

def format_number(num, prec=3):
    '''Returns the formated number up to specific precision'''
    fmt = '{:,.' + str(prec) + 'f}'
    s = fmt.format(num)
    return s.rstrip('0').rstrip('.')

def fisher_exact(table):
    table = np.array(table)
    r.assign('cont', np2r(table))
    if np.max(table) < 1000 or (table.shape[0] == 2 and table.shape[1] == 2):
        r('res <- fisher.test(cont)')
    else:
        r('res <- fisher.test(cont, simulate.p.value=TRUE)')
    return float(r('res')[0][0])

def table_one(df, grpvar=None, catvars=None):
    # convert boolean type to int --> remains are float or str
    df.replace({False: 0, True: 1}, inplace=True)

    if grpvar is not None:
        grp_names = df[grpvar].unique()

    # for csv file
    rows = []

    # Generate table header
    tabs = ['', 'Total']
    if grpvar is not None:
        for grp_name in grp_names:
            tabs.append(f'{grpvar}={grp_name}')
        if len(grp_names) > 1:
            tabs.append('P-value')
            tabs.append('Test')
    rows.append(tabs)

    tabs = ['n']
    if grpvar is not None:
        tabs.append(format_number(sum(~df[grpvar].isnull())))
        for grp_name in grp_names:
            tabs.append(format_number(sum(df[grpvar] == grp_name)) + ' (' + format_number(np.mean(df[grpvar] == grp_name) * 100, 1) + '%)')
    else:
        tabs.append(format_number(len(df)))
    rows.append(tabs)

    # Generate statistics for each variable
    for col in df.columns:
        if col == grpvar:
            continue
        try:
            pd.to_numeric(df[col])
            isstr = False
        except:
            isstr = True

        unique_values = sorted(df.loc[~df[col].isnull(), col].unique())  # unique values

        iscat = len(unique_values) < 8
        if catvars:
            if col in catvars:
                iscat = True

        if isstr and not iscat:
            continue

        if iscat:  # categorical variables --> represents as count (percent)
            if grpvar is not None:  # create cross table (value x grp)
                xtab = pd.crosstab(df[col], df[grpvar]).fillna(0)
                pval = None
                if len(grp_names) > 1 and min(xtab.shape) >= 2: # NEJM requires Exact method for all categorical variables
                    if (xtab > 5).all(axis=None):  # if there is an incidence < 5
                        pval = stat.chi2_contingency(xtab)[1]
                        test_name = 'Chi-square'
                    else:
                        pval = fisher_exact(xtab.T.values)
                        test_name = 'Fisher\'s exact'

            is_binary = (len(unique_values) == 2) and (unique_values[0] == 0 and unique_values[1] == 1)
            if is_binary:  # binary
                # print total
                tabs = [col, format_number(sum(df[col] == 1)) + ' (' + format_number(np.mean(df[col] == 1) * 100, 1) + '%)']
                if grpvar is not None: # print group values
                    for grp_name in grp_names:
                        grp_mask = (df[grpvar] == grp_name)
                        tabs.append(format_number(sum(df.loc[grp_mask, col])) + ' (' + format_number(np.mean(df.loc[grp_mask, col])*100, 1) + '%)')
                    if pval is not None:
                        tabs.append(format_pval(pval))
                        tabs.append(test_name)
                rows.append(tabs)
            else:
                for uval in unique_values:
                    # print total
                    tabs = [f'{col}={uval}', format_number(sum(df[col] == uval)) + ' (' + format_number(np.mean(df[col] == uval) * 100, 1) + '%)']
                    if grpvar is not None: # print group values
                        for grp_name in grp_names:
                            grp_mask = (df[grpvar] == grp_name)
                            tabs.append(format_number(sum(df.loc[grp_mask, col] == uval)) + ' (' + format_number(np.mean(df.loc[grp_mask, col] == uval) * 100, 1) + '%)')
                        if pval is not None:
                            if uval == unique_values[0]:
                                tabs.append(format_pval(pval))
                                tabs.append(test_name)
                    rows.append(tabs)

        else:  # continuous variables --> represents as mean (SD)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                isnorm = stat.shapiro(df[col])[1] > 0.05  # check if it is normal distribution
            if isnorm:  # normal distribution
                # print total
                m = df[col].mean()
                s = df[col].std()
                tabs = [col, f'{m:.3f} ({s:.3f})']

                if grpvar is not None:
                    # extract group values
                    grp_vals = []
                    for grp_name in grp_names:
                        a = df.loc[df[grpvar] == grp_name, col]
                        grp_vals.append(a[~a.isnull()])

                    # print group values
                    for igrp in range(len(grp_vals)):
                        m = grp_vals[igrp].mean()
                        s = grp_vals[igrp].std()
                        tabs.append(f'{m:.3f} ({s:.3f})')

                    # print stats
                    if len(grp_names) == 2:
                        equal_var = stat.levene(grp_vals[0], grp_vals[1])[1] > 0.05  # levene
                        pval = stat.ttest_ind(grp_vals[0], grp_vals[1], equal_var=equal_var)[1]
                        test_name = 'T-test'
                    else:  # 3 or more groups -> anova
                        equal_var = stat.levene(*grp_vals)[1] > 0.05  # levene + homoscedasticity
                        if equal_var:
                            pval = stat.f_oneway(*grp_vals)[1]
                            test_name = 'One-way ANOVA'
                        else:
                            pval = stat.kruskal(*grp_vals)[1]
                            test_name = 'Kruskal-Wallis'
                        tabs.append(format_pval(pval))
                        tabs.append(test_name)
            else:  # non-normal
                # print total
                m = df[col].median()
                q1 = df[col].quantile(0.25)
                q2 = df[col].quantile(0.75)
                tabs = [col, format_number(m, 3) + ' (' + format_number(q1, 3) + '-' + format_number(q2, 3) + ')']

                if grpvar is not None:
                    # extract group values
                    grp_vals = []
                    for grp_name in grp_names:
                        a = df.loc[df[grpvar] == grp_name, col]
                        grp_vals.append(a[~a.isnull()])

                    # print group value
                    for igrp in range(len(grp_vals)):
                        m = grp_vals[igrp].median()
                        q1 = grp_vals[igrp].quantile(0.25)
                        q2 = grp_vals[igrp].quantile(0.75)
                        tabs.append(format_number(m, 3) + ' (' + format_number(q1, 3) + '-' + format_number(q2, 3) + ')')

                    # print stats
                    if len(grp_vals) == 2:
                        pval = stat.mannwhitneyu(grp_vals[0], grp_vals[1], alternative='two-sided')[1]
                        test_name = 'Mann-Whitney'
                    elif len(grp_vals) > 2:  # > 3 groups
                        pval = stat.kruskal(*grp_vals)[1]
                        test_name = 'Kruskal-Wallis'
                    tabs.append(format_pval(pval))
                    tabs.append(test_name)

            rows.append(tabs)

    return pd.DataFrame(rows)

# read data
df = pd.read_csv('https://api.vitaldb.net/cases')

# add columns
df['opdur'] = df['opend'] - df['opstart']
df['anedur'] = df['aneend'] - df['anestart']
df['hospdur'] = df['dis'] - df['adm']

# remove columns
df.drop(columns=['opstart', 'opend', 'anestart', 'aneend', 'dis', 'adm', 'caseid'], inplace=True)
df = df.loc[:, ~df.columns.str.endswith('id')]

# create table one
df_results = table_one(df, 'department')
#df_results = table_one(df, 'death_inhosp')

# save and print results
df_results.to_csv('table1.csv', index=False, header=False)
print(df_results)