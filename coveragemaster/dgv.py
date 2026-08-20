"""Database of Genomic Variants overlap queries."""
import os
import warnings

import numpy as np
import pandas as pd


class DGVExplorer:
    """Overlap and frequency lookup against a DGV tabular file."""

    def __init__(self, file_path, filt_cols, chr_idx_path):
        self.file_path = file_path
        self.filt_cols = filt_cols
        self.chr_idx_path = chr_idx_path
        self.data = self._read_csv()
        if not os.path.exists(self.chr_idx_path):
            self.chr_indices = self._create_chr_indices()
        else:
            self.chr_indices = np.load(self.chr_idx_path, allow_pickle=True).item()

    def _read_csv(self):
        return pd.read_csv(self.file_path, delimiter="\t", usecols=self.filt_cols, header=0)

    def _create_chr_indices(self):
        chr_indices = {}
        for chrom in self.data["#chrom"].unique():
            indices = self.data.index[self.data["#chrom"] == chrom]
            chr_indices[chrom] = (indices[0], indices[-1])
        np.save(self.chr_idx_path, chr_indices, allow_pickle=True)
        return chr_indices

    def get_overlap(self, chrom, alpha_0, beta_0):
        overlaps = self.data.iloc[
            self.chr_indices[chrom][0]: self.chr_indices[chrom][1] + 1
        ]
        idx = overlaps["chromStart"].searchsorted(beta_0, side="right").item()

        if idx > 0:
            overlaps = overlaps.iloc[:idx]
        elif overlaps["chromStart"].iloc[0] <= beta_0:
            overlaps = overlaps.iloc[[0]]
        else:
            warnings.warn(f"No overlap found [chrom:{chrom}, a:{alpha_0}, b:{beta_0}]")
            return None

        if len(overlaps) > 1:
            idx2 = overlaps["chromStart"].searchsorted(alpha_0, side="left").item()
            if idx2 != len(overlaps):
                data_in = overlaps.iloc[idx2:]
                data_out = overlaps.iloc[:idx2]
                overlaps = pd.concat([
                    data_out[
                        (data_out["chromStart"] < alpha_0) & (data_out["chromEnd"] >= alpha_0)
                    ],
                    data_in,
                ])
            else:
                overlaps = overlaps[
                    (overlaps["chromStart"] < alpha_0) & (overlaps["chromEnd"] >= alpha_0)
                ]

        if len(overlaps) == 0:
            warnings.warn(f"No overlap found [chrom:{chrom}, a:{alpha_0}, b:{beta_0}]")
            return None
        return overlaps

    @staticmethod
    def get_freq(data):
        return {
            "freq_g": data["observedGains"] / data["sampleSize"],
            "freq_l": data["observedLosses"] / data["sampleSize"],
            "name": data["name"],
        }
