from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import scipy
from module.getters import (
    getCompoundAndRatiosDf,
    getAggregateStatsDf,
    getExperimentalInfo,
    getTreatmentMapping,
    getRegionSubclassification,
    updateQuantitativeStats,
    getQuantitativeStats,
)
from module.utils import (
    askMultipleChoice,
    askSelectParameter,
    askYesorNo,
    figure_cache,
    subselectDf,
)

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt


def pca(
        filename,
        experiment=None,
        compounds=None,
        regions=None,
        from_scratch=None,
        threshold_vairance = 0.1 #explained vairance of PC > treshold vairance
    ):
    data = getCompoundAndRatiosDf(filename)
    experimental_info = getExperimentalInfo(filename)
    exit_loop = False
    while not exit_loop:
        # We use keyword params here even though they actually are mandatory for the decorator
        singlePCA(
            filename,
            experiment=experiment
            or askMultipleChoice("Select experiment", experimental_info.keys()),
            compounds=compounds or askSelectParameter(data, "compound"), #currently taking list only #REMI improve/make same as others
            regions=regions or askSelectParameter(data, "region"),
            from_scratch=from_scratch,
        )
        exit_loop = askYesorNo("Exit?")


@figure_cache("pca")
def singlePCA(
        filename,
        experiment=None,
        compounds=None,
        regions=None,
        from_scratch=None,
        threshold_vairance=0.1
):
    #fetch df subselection and restructre df
    df = subselectDf(getCompoundAndRatiosDf(filename), {"experiment": experiment, "compound":compounds, "region":regions})
    df['compound_region'] = df['compound'] + '_' + df['region']
    df_X = df.pivot( index = [ 'mouse_id', 'treatment', 'color'], columns=['compound_region'], values='value')
    df_X.reset_index(inplace=True)
    df_X.columns.name = None


    #handeling and reporting NaN values
    df_X = df_X.replace(0, np.nan) #FIX ME this should not be nessary but the df.pkl should have nan not 0
    nan_rows = []
    for index, row in df_X.iterrows():
        if row.isna().any():
            columns_with_nan = row.index[row.isna()] # Get the column names with NaN values
            nan_rows.append((index, columns_with_nan))
    for row_index, cols in nan_rows: # Print the rows with NaN values and the corresponding column names
        mouse=df_X.loc[row_index,'mouse_id']
        print(f"Mouse id {mouse} has NaN values in columns: {', '.join(cols)}")
        df_X=df_X.drop(row_index)
        print(f"... removed mouse {mouse}")


    mouse_id = df_X[['mouse_id', 'treatment', 'color']] #  non-value columns
    features = df_X.drop(['mouse_id', 'treatment', 'color'], axis=1) # values with columns=features i.e. DA_OF

    # Standardize the features
    scaler = StandardScaler()
    features_standardized = scaler.fit_transform(features)

    # Perform PCA
    pca = PCA()
    pca_result = pca.fit_transform(features_standardized)


    # create df with pca results, mouse_id df back
    df_pca = pd.DataFrame(data=pca_result, columns=[f'PC{i+1}' for i in range(features_standardized.shape[1])])
    df_pca = pd.concat([mouse_id, df_pca], axis=1)

    #PC's>0.1 explained vairance 
    explained_variance = pca.explained_variance_ratio_
    # Filter PCs with explained variance > 0.1
    selected_pcs = [f'PC{i+1}' for i, variance in enumerate(explained_variance) if variance > threshold_vairance]


    # Create a figure 
    fig_PCs, ax_PCs = plt.subplots(figsize=(8, 6))

    # Extract the explained variances for the selected PCs
    selected_variances = [explained_variance[int(pc[2:]) - 1] for pc in selected_pcs]

    # Plot the PCs with explained variance for each PC > 0.1
    ax_PCs.bar(range(1, len(selected_pcs) + 1), selected_variances, align='center')
    # ax_PCs.xticks(range(1, len(selected_pcs) + 1), selected_pcs)
    ax_PCs.set_xticks(range(1, len(selected_pcs) + 1))
    ax_PCs.set_xlabel('Principal Component')
    ax_PCs.set_ylabel('Explained Variance')
    ax_PCs.set_title('Principal Components with Explained Vairance > 0.1')
    plt.show() #REMI as i cant save multiple figures i just have it spat out in the notebook seems suboptimal

    # Iterate over the selected PCs print LOADINGS
    for pc_label in selected_pcs:
        pc_index = int(pc_label[2:]) - 1  # Extract the PC index (e.g., PC1 -> 0)
        loadings = pca.components_[pc_index]
        loadings_df = pd.DataFrame({'Factor': features.columns, 'Loading': loadings})

        # Sort the factors by the absolute magnitude of loading in descending order
        loadings_df['Absolute Magnitude'] = loadings_df['Loading'].abs()
        loadings_df = loadings_df.sort_values(by='Absolute Magnitude', ascending=False)
        loadings_df = loadings_df.drop(columns=['Absolute Magnitude'])

        # Print the loadings table for the current PC
        print(f"Loadings (Factors) for {pc_label} (Ordered by Absolute Magnitude):")
        print(loadings_df)

    #PLOT PCA 
    fig_pca, ax_pca = plt.subplots(figsize=(10, 6))
    # Calculate the sum of variances for PC1 and PC2
    summed_variance = selected_variances[0] + selected_variances[1]

    # colors based on color column
    color_dict = df_pca.set_index('treatment')['color'].to_dict()
    for treatment, color in df_pca.groupby('treatment'):
        ax_pca.scatter(
            color['PC1'],
            color['PC2'],
            label=treatment,
            color=color_dict[treatment],
            edgecolor='black',  # Add a black outline
            alpha=0.5
        )

    ax_pca.set_xlabel('Principal Component 1')
    ax_pca.set_ylabel('Principal Component 2')
    ax_pca.legend()
    ax_pca.set_title(f'PCA (Summed Variance: {summed_variance:.2f})')

    
    return  fig_pca 
