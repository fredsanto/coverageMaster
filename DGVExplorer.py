import pandas as pd
import numpy as np

import os
import warnings


class DGVExplorer:

    def __init__(self, file_path, filt_cols, chr_idx_path):
        """
        Initializes a DVGExplorer object for reading and manipulating dgv files

        :param file_path: Path to the .csv file
        :param filt_cols: Columns to select/filter in the dataframe
        :param chr_idx_path: Path to the .npy file containing the chromosome indices

        :type file_path: str
        :type filt_cols: list
        :type chr_idx_path: str
        """

        self.file_path = file_path
        self.filt_cols = filt_cols
        self.chr_idx_path = chr_idx_path

        self.data = self._read_csv()

        if not os.path.exists(self.chr_idx_path):
            self.chr_indices = self._create_chr_indices()
        else:
            # TODO: check for encoding='latin1' when saving in python2 but loading in python3
            self.chr_indices = np.load(self.chr_idx_path, allow_pickle=True).item()

        # self.overlaps = None

    def _read_csv(self):
        """
        Reads the DGV .csv file
        :return: A Dataframe containing the DGV data
        """
        return pd.read_csv(self.file_path, delimiter='\t', usecols=self.filt_cols, header=0)

    def _create_chr_indices(self):
        """
        Creates and stores a dictionary with the (start, end) line indices of each chromosome inside the Dataframe

        :return: Dictionary containing the (start, end) line indices  of each chromosome inside the Dataframe
        :rtype: dict
        """
        chr_in_data = self.data["#chrom"].unique()
        chr_indices = {}

        for chrom in chr_in_data:
            indices = self.data.index[self.data["#chrom"] == chrom]
            chr_indices[chrom] = (indices[0], indices[-1])

        np.save(self.chr_idx_path, chr_indices, allow_pickle=True)

        return chr_indices

    def get_overlap(self, chrom, alpha_0, beta_0):
        """
        Finds the overlap of CM call (a0,b0) with the DGV calls (a,b).
        1st. keep overlaps if a<=b0 (search sorted)
        2nd. from 1st case, keep overlaps if a0<=a (<=b) (search sorted)
        3rd. keep overlaps if a0>a and b>=b0 (simple indexing)

        :param chrom:
        :param alpha_0:
        :param beta_0:
        :return: The overlapping data (if exist)
        :rtype: pd.DataFrame
        """

        overlaps = self.data.iloc[self.chr_indices[chrom][0]:self.chr_indices[chrom][1] + 1]

        # log search to get data indices for which b0 >= a to exclude definite non-overlaps
        alpha_smaller_equal_beta_0_idx = overlaps["chromStart"].searchsorted(beta_0, side='right').item()

        # Reduce the data and ensure that we do not have a degenerate case
        if alpha_smaller_equal_beta_0_idx > 0:
            overlaps = overlaps.iloc[0:alpha_smaller_equal_beta_0_idx]
        # if index = 0, then the 1st line of DGV file has a probable overlapping call. If DGV(start) <= bo, then it s the only overlap.
        elif overlaps["chromStart"][alpha_smaller_equal_beta_0_idx] <= beta_0:
            overlaps = overlaps.iloc[0]
        else:
            # if index = 0, then there is no overlap
            warnings.warn("No overlap found [chrom:%s, a:%d, b:%d]" % (chrom, alpha_0, beta_0))
            return None

        # In the degenerate case, we will only have a single overlap hence we do not need to check any other interval.
        # Otherwise, we should get the rest of the intervals
        if len(overlaps) > 1:
            # log search to get data indices for which a0 >= a
            alpha_greater_equal_alpha_0_idx = overlaps["chromStart"].searchsorted(alpha_0, side='left').item()

            # Ignore degenerate case where no gene starts within [a0, b0]
            # if index != last line of file:
            if alpha_greater_equal_alpha_0_idx != len(overlaps):

                data_alpha_in_a0_b0 = overlaps.iloc[alpha_greater_equal_alpha_0_idx:len(overlaps)]
                data_alpha_smaller_than_a0 = overlaps.iloc[0:alpha_greater_equal_alpha_0_idx]

                overlaps = pd.concat([
                    data_alpha_smaller_than_a0[(data_alpha_smaller_than_a0["chromStart"] < alpha_0) & (
                                data_alpha_smaller_than_a0["chromEnd"] >= alpha_0)],
                    data_alpha_in_a0_b0
                ])
            else:
                overlaps = overlaps[(overlaps["chromStart"] < alpha_0) & (overlaps["chromEnd"] >= alpha_0)]

        if len(overlaps) == 0:
            warnings.warn("No overlap found [chrom:%s, a:%d, b:%d]" % (chrom, alpha_0, beta_0))
            return None

        # self.overlaps = overlaps
        return overlaps

    @staticmethod  # static coz we do not use anything related to self
    def get_freq(data):
        freq_g = data['observedGains'] / data['sampleSize']
        freq_l = data['observedLosses'] / data['sampleSize']
        name = data['name']

        # print(["{0:.3f}".format(i) for i in freq_g.tolist()])
        # print(["{0:.3f}".format(i) for i in freq_l.tolist()])

        return {"freq_g": freq_g,
                "freq_l": freq_l,
                "name": name
                }

#
# filt_columns = ["#chrom",
#                 "chromStart",
#                 "chromEnd",
#                 "name",
#                 "varType",
#                 "supportingVariants",
#                 "sampleSize",
#                 "observedGains",
#                 "observedLosses"]
#
# DGV_file_path = "/users/mrapti/scratch/melivoia/DGV/DGVtable_hg19"
# chr_idx_path = "/users/mrapti/scratch/melivoia/DGV/DGVtable_hg19_chr_indices.npy"
#
# explorer = DGVExplorer(file_path=DGV_file_path,
#                        filt_cols=filt_columns,
#                        chr_idx_path=chr_idx_path)
'''
overlaps = explorer.get_overlap(chrom="chr1",
                                alpha_0=1,
                                beta_0=1000)
print (overlaps)

if overlaps is not None:
    freq_results = explorer.get_freq(overlaps)
    print (freq_results)

    print (overlaps['varType'])
'''