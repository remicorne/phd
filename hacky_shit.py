df = getCompoundDf(filename)
def passes_pvalue_threshold(threshold):
    return lambda x, y: pearsonr_pval(x, y) >= threshold

def pearsonr_pval(x,y):
        return scipy.stats.pearsonr(x,y)[1]
        
def pearsonr_r(x,y):
        return scipy.stats.pearsonr(x,y)[0]


BR_column_order = ['OF', 'PL', 'CC', 'IC', 'M', 'SJ', 'SL1', 'SR1', 'SL6', 'SR6', 'AC', 'V' , 
            'Am', 'dH', 'vH', 'NAc', 'VM', 'DM', 'VL', 'DL', 'MD', 'VPL', 'VPR', 'DG', 
            'Y', 'SC', 'SN', 'VTA', 'DR', 'MR', 'CE' ] 


#within a compound between BRs
for treatment_compound, treatment_compound_groupby_df in df.groupby(by =['treatment', 'compound']):
    treatment, compound = treatment_compound
    #create df_to_correlate: index as mouse number, column as BR, ng_mg and NaN where data missing 
    # print(f'reordered BR_list for {treatment} within {compound} : {BR_list}')

    df_to_correlate = treatment_compound_groupby_df.pivot_table(values='ng_mg', index=df[['mouse_id']], columns='BR').filter(BR_column_order)
    corr_matrix_BR = df_to_correlate.corr('pearson', min_periods=5).dropna(axis=0, how='all').dropna(axis=1, how='all')
    p_value_mask = df_to_correlate.corr(method=passes_pvalue_threshold(0.05), min_periods=5).dropna(axis=0, how='all').dropna(axis=1, how='all').astype(bool)

    #plot shit #this runs but the figures dont get spat out? REMI WHYYYY
    fig, ax = plt.subplots(figsize=(16, 10))
    # mask = np.invert(np.tril(p_value_mask < 0.05)) #generalise p_value = 0.05

    heatmap = sns.heatmap(corr_matrix_BR, vmin=-1, vmax=1,
                            annot=True, cmap='BrBG', mask=p_value_mask, annot_kws={"size": 8})
    heatmap.set_title(
        compound+'  in  ' + treatment, fontdict={'fontsize': 18}, pad=12)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        horizontalalignment='right')
    ax.set_ylabel('')
    ax.set_xlabel('')
    fig.savefig('test1.png')
    plt.close('all')