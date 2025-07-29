# ----------------------------------------------------------------------------
# Copyright (c) 2025, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os

import pandas as pd
import qiime2

from ._format import GUNCResultsDirectoryFormat


def gunc_results_directory_format_to_metadata(
    fmt: GUNCResultsDirectoryFormat,
) -> qiime2.Metadata:
    """Transform GUNCResultsDirectoryFormat to qiime2.Metadata.
    
    Parameters
    ----------
    fmt : GUNCResultsDirectoryFormat
        GUNC results directory format to transform.
        
    Returns
    -------
    qiime2.Metadata
        Metadata object with MAG IDs as index and optional sample_id column.
    """
    file_dict = fmt.file_dict()
    
    # Collect all dataframes from maxCSS files
    dataframes = []
    
    for sample_id, directory in file_dict.items():
        # Find the maxCSS file in this directory
        pattern = os.path.join(directory, "GUNC.*.maxCSS_level.tsv")
        maxcss_files = glob.glob(pattern)
        
        if not maxcss_files:
            # Skip directories without maxCSS files
            continue
            
        # Use the first (and should be only) maxCSS file
        maxcss_file = maxcss_files[0]
        
        # Read the TSV file
        df = pd.read_csv(maxcss_file, sep='\t')
        
        # If we have sample data (partitioned), add sample_id column
        if sample_id:  # sample_id is not empty string
            df['sample_id'] = sample_id
            
        dataframes.append(df)
    
    # Concatenate all dataframes
    if not dataframes:
        raise ValueError("No GUNC maxCSS results found in the directory format")
        
    combined_df = pd.concat(dataframes, ignore_index=True)
    
    # Set MAG ID as index with name "id"
    combined_df.set_index('genome', inplace=True)
    combined_df.index.name = 'id'
    
    # Create and return Metadata object
    return qiime2.Metadata(combined_df)