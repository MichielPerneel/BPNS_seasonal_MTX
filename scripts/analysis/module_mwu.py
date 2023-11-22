import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def qvals_bh(pvals: np.array):
    return multipletests(pvals, method="fdr_bh")[1]


def enrich_mwu(ordered_list, annotation: pd.DataFrame,
               min_size=5, alternative="two-sided", print_progress=0):
    """
    Check if a KEGG pathway occurs at the top/bottom of an ordered list of KEGG KOs.

    :param ordered_KO_list: List of KEGG KO ids associated with a module, ordered by some metric (e.g. moduleMembership).
    :param annotation: columns: "KEGG_ko", "annotation"
    :param min_size: the minimum amount of KOs in annotation category (that are also in ordered list)
        that we need in order to actually run the MWU test.
        I. e. smaller gene sets are skipped.
    :param alternative: passed to `mannwhitneyu`.
    :param print_progress: print a progress message every 'n' tests.
    :return: A data frame with columns:
        - annotation: the annotation category that was tested
        - p_val: the two-sided p-value of the MWU test:
          the more a GO term is located at the top/bottom, the lower the p-value.
        - q_val: the BH-adjusted p-value.
        - U: the MWU test-statistic.
        - size: the size of the overlap between the number of KOs belonging to a annotation category and all KO IDs in the list
        - estimate: the probability that a random KO of the given category
        is ranked higher than a random KO of the ordered list.
        - direction: can take three values: "top"/"mid"/"bot", if the estimate is bigger/equal/smaller than 0.5.
          For example, "top" means a KO id is located closer to the to of the given list.
    """
    # Get the indices of all identifiers in the ordered list.
    all_indices = list(range(len(ordered_list)))

    # Group the KEGG KOs by their associated pathway.
    groups = annotation.groupby("annotation")
    n_groups = len(groups)
    results = []
    
    # Loop over each pathway.
    for i, (annotation, group) in enumerate(groups):
        # Print a progress message if print_progress is set and i is a multiple of print_progress
        if print_progress and i % print_progress == 0:
            print("{:6d} / {:6d} ({:05.2f}%)".format(i, n_groups, i/n_groups * 100))

        # Get the KO ids of the current pathway.
        pathway_KO_IDs = group.KEGG_ID.values
        
        # Find the indices of the KO ids in the ordered list that are also in the current pathway.
        indices_of_annotation = np.flatnonzero(np.in1d(ordered_list, pathway_KO_IDs))

        # Only perform the MWU test if the number of overlapping genes is at least min_size
        if len(indices_of_annotation) >= min_size:
            # Perform the Mann-Whitney U test between the indices of all genes in the ordered list and
            # the indices of the genes in the current GO category
            u2, p_val = mannwhitneyu(x=all_indices, y=indices_of_annotation, alternative=alternative)
            results.append({"annotation": annotation, "p_val": p_val, "U": u2, "size": len(indices_of_annotation)})

    # Record the results in a dictionary
    df = pd.DataFrame.from_records(results)
    # Add the corrected p-values
    df.insert(column="q_val", value=qvals_bh(df["p_val"]), loc=2)
    df["estimate"] = df["U"]/(df["size"] * len(all_indices))
    df["direction"] = df["estimate"].apply(lambda e: "mid" if e == 0.5 else ("top" if e > 0.5 else "bot"))
    df["_extremity"] = 0.5 - abs(0.5 - df["estimate"])
    df.sort_values(by=["p_val", "_extremity"], inplace=True)
    df.drop(columns="_extremity", inplace=True)
    return df


def main():
    ordered_list = ['a', 'c', 'b', 'f', 'e', 'd']
    categs = pd.DataFrame(dict(gene_id=['b', 'a', 'f'], annotation=['x', 'x', 'x']))
    res = enrich_mwu(ordered_list, categs, min_size=0)
    print(res)


if __name__ == '__main__':
    main()