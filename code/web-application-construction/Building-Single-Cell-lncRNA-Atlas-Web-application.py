#!/usr/bin/env python
# coding: utf-8


# In[ ]:


"""
Python script to build Single-Cell lncRNA Atlas Web application
Hyunsu An @ Gwangju Institute of Science and Technology
# last modified 2024-02-29 10:51:17 
"""
'''python version == Python 3.9.13 '''
# pip install biobookshelf==0.2.10
''' 
initialize 
'''
from biobookshelf.main import *
from biobookshelf import *
import biobookshelf as bk

path_folder_wd = '/home/project/Single_Cell_lncRNA_Database/data/'
def Read_Expression_Data( path_file, dtype = 'float16' ) :
    # store float values as 'float16' to reduce memory footprint
    dict_col_to_dtype = { "": 'str' }
    with open( path_file ) as file :
        l_col = file.readline( ).strip( ).split( '\t' )
    for col in l_col :
        dict_col_to_dtype[ col ] = dtype
    return pd.read_csv( path_file, sep = '\t', index_col = 0, dtype = dict_col_to_dtype )


# In[ ]:


''' 
load expression data 
# 2023-06-28 15:50:41 
'''

"""
prepare gene expression matrix
"""
df_expression_mTissues6 = pd.read_csv( f'{path_folder_wd}processed_data/atlas/multi_tissue/whole_tissue_exprssion_matrix.csv', index_col = 0 ) # read data
df_expression_aging_kidney = pd.read_csv( f'{path_folder_wd}processed_data/20240226_update/atlas_R_object/aging_kidney_metadata/aging_kidney_matrix.csv', index_col = 0 ) 
df_expression_mTissues6.index.name = 'name_gene'
df_expression_aging_kidney.index.name = 'name_gene'
# merge the two dense matrix 
# reset index
df_expression_mTissues6.reset_index( inplace = True )
df_expression_aging_kidney.reset_index( inplace = True )
df_expression = pd.merge( df_expression_mTissues6, df_expression_aging_kidney, on = [ "name_gene" ], how = 'outer' ) # merge to create a new data
df_expression.fillna( 0, inplace = True ) # fill na values
df_expression.set_index( 'name_gene', inplace = True ) # set index

"""
pre-process gene expression matrix and perform sanity checks
"""
df_expression.columns = list( e.replace( '.', '-' ) for e in df_expression.columns.values ) # change cell names for compatibility
df_expression.sort_index( inplace = True ) # sort by gene_names to increase a chance to reuse the fetched data (clustering similar genes) # sort by gene names to increase the probability of finding similar genes within a chunk of file
print( "check whether id_cells of expression dataset are unique... ", len( LIST_COUNT( df_expression.columns.values ) ) == 0 )
print( "check whether the expression data contains nan values... ", pd.isnull( df_expression ).sum( ).sum( ) != 0 )
s_expr_mean = df_expression.mean( axis = 1 ).sort_values( )
s_expr_max = df_expression.max( axis = 1 ).sort_values( )


# In[ ]:


'''
retrieve UMAP coordinate data and basic cell quality information of each tissue 
# modified 2024-02-28 15:09:36 
'''
# previous annotations
df_cells = pd.read_csv( path_folder_wd + 'processed_data/atlas/multi_tissue/whole_tissue_celltype_coordination.csv', index_col = 0 ) # load global UMAP coordinates
df_cells.rename( columns = { 'UMAP1' : "UMAP_1", 'UMAP2' : 'UMAP_2', 'ident' : 'cell_type', 'orig.ident' : 'tissue', 'percent.mt' : 'percent_mt' }, inplace = True )
df_cells_tissue = pd.concat( list( pd.read_csv( path_file, index_col = 0 ) for str_tissue, path_file in bk.PD_Select( GLOB_Retrive_Strings_in_Wildcards( f"{path_folder_wd}processed_data/atlas/multi_tissue/*_celltype_coordination.csv" )[ [ 'wildcard_0', 'path' ] ], wildcard_0 = 'whole_tissue', deselect = True ).values ) )    
df_cells_tissue = df_cells_tissue[ [ 'ident', 'UMAP1', 'UMAP2', ] ]
df_cells_tissue.rename( columns = { 'UMAP1' : "tissue_UMAP_1", 'UMAP2' : 'tissue_UMAP_2', 'ident' : 'tissue_cell_type', }, inplace = True )
df_cells = df_cells.join( df_cells_tissue ) # combine global and tissue level data

df_cells.tissue_cell_type = df_cells.tissue_cell_type.apply( MAP.Map( { "Naïve_Cd4" : "Naive_Cd4", 'Naïve_Cd8' : "Naive_Cd8", 'Naïve T' : 'Naive T', 'αβ T' : 'alpha beta T', 'CD8 αα' : 'CD8 alpha alpha' } ).a2b_if_mapping_available_else_Map_a2a ) # adjust annotations # extended ASCII character is currently unsupported

# add aging kidney data
df_cells_aging_kidney = pd.read_csv( path_folder_wd + 'processed_data/20240226_update/atlas_R_object/aging_kidney_metadata/aging_kidney_metadta.csv', index_col = 0 ) # read new data
df_cells_aging_kidney.rename( columns = { 'UMAP_1' : "tissue_UMAP_1", 'UMAP_2' : 'tissue_UMAP_2', 'ident' : 'tissue_cell_type', 'percent.mt' : 'percent_mt', 'information' : 'orig.ident' }, inplace = True )

# add columns for integrations
# create new columns
df_cells[ 'origin' ] = df_cells[ 'tissue' ]
df_cells[ 'age' ] = '6-8 weeks (young mice)'
df_cells_aging_kidney[ 'origin' ] = df_cells_aging_kidney[ 'tissue' ]

for name_col in [ 'cell_type', 'UMAP_1', 'UMAP_2' ] :
    df_cells_aging_kidney[ name_col ] = df_cells_aging_kidney[ f"tissue_{name_col}" ]
    
# modify tissue column    

l = [ ]
for origin, age in df_cells_aging_kidney[ [ 'tissue', 'age' ] ].values :
    l.append( f"{origin.capitalize( )}-enriched ({age} mice)" )
df_cells_aging_kidney[ 'tissue' ] = l

# create consistent labels for 'age' column
dict_map = { 'old' : '84 weeks (old mice)', 'young' : '8 weeks (young mice)' }
df_cells_aging_kidney[ 'age' ] = list( dict_map[ e ] for e in df_cells_aging_kidney[ 'age' ].values ) 

# compose 'df_cells'
df_cells[ 'tissue' ] = list( e.capitalize( ) for e in df_cells.tissue.values )
df_cells = pd.concat( [ df_cells, df_cells_aging_kidney[ df_cells.columns.values ] ] ) # compose 'df_cells'
df_cells[ 'origin' ] = list( e.capitalize( ) for e in df_cells.origin.values )


# In[ ]:


def Round_Float( df, l_col_scientific_notations, l_col_typical_notation, n_significant_digits_scientific_notation = 3, n_significant_digits_typical_notation = 3, inplace = False ) :
    ''' round float with a given number of significant digits and convert floating point numbers to strings '''
    if not inplace :
        df = deepcopy( df )
    str_format_scientific_notation = "{:." + str( int( n_significant_digits_scientific_notation ) ) + "e}"
    str_format_typical_notation = '{:.' + str( int( n_significant_digits_typical_notation ) ) + 'f}'
    for col in l_col_scientific_notations :
        df[ col ] = list( '' if np.isnan( value ) else str_format_scientific_notation.format( value ) for value in df[ col ].values )
    for col in l_col_typical_notation :
        df[ col ] = list( '' if np.isnan( value ) else str_format_typical_notation.format( value ) for value in df[ col ].values )
    return df


# In[ ]:


"""
retrieve cellmarkers from 'FindAllMarkers' program 
# 2024-02-28 15:12:45 
"""
l_df = [ ]
# mTissues6 dataset
for name_tissue in [ 'heart', 'intestine', 'kidney', 'lung', 'liver', 'thymus' ] :
    df = pd.read_csv( f"{path_folder_wd}processed_data/atlas/marker/{name_tissue}_marker.csv", index_col = 0 )
    df.columns = [ 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster', 'gene' ]
    df[ 'tissue' ] = name_tissue.capitalize( )
    l_df.append( df )
# aging kidney data
df = pd.read_csv( f"{path_folder_wd}processed_data/atlas/marker/aging_kidney_marker.csv", index_col = 0 )
df.columns = [ 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster', 'gene' ]
for name_tissue in [ 'Glom-enriched (old mice)', 'Tubule-enriched (old mice)', 'Glom-enriched (young mice)', 'Tubule-enriched (young mice)' ] :
    df = df.__deepcopy__( ) # copy the dataframe
    df[ 'tissue' ] = name_tissue # use the same cell type for all aging kidney data
    l_df.append( df )
    
df_cell_marker = pd.concat( l_df, ignore_index = True ) # reset index
df_cell_marker = Round_Float( df_cell_marker, l_col_scientific_notations = [ 'p_val', 'p_val_adj' ], l_col_typical_notation = [ 'avg_log2FC', 'pct.1', 'pct.2' ], n_significant_digits_scientific_notation = 3, n_significant_digits_typical_notation = 3 )
df_cell_marker.cluster = df_cell_marker.cluster.apply( MAP.Map( { "Naïve_Cd4" : "Naive_Cd4", 'Naïve_Cd8' : "Naive_Cd8", 'Naïve T' : 'Naive T', 'αβ T' : 'alpha beta T', 'CD8 αα' : 'CD8 alpha alpha' } ).a2b_if_mapping_available_else_Map_a2a ) # adjust annotations to match that of 'df_cells.tissue_cell_type' # 2023-06-30 17:24:45 
# df_cell_marker.cluster = df_cell_marker.cluster.apply( MAP.Map( { 'rd T' : "gd T", 'Cd8 aa' : 'CD8 aa' } ).a2b_if_mapping_available_else_Map_a2a ) # adjust annotations to match that of 'df_cells.tissue_cell_type'

# index cellmarker data
WEB.Index_and_Base64_Encode( df_to_be_indexed = df_cell_marker, l_col_index = [ 'tissue', 'cluster' ], path_prefix_output = f"{path_folder_wd}processed_data/web_application/base64/cell_marker", path_folder_temp = '/tmp/' )


# In[ ]:


''' export expression data gene by gene '''
# shutil.rmtree( f"{path_folder_wd}processed_data/expression/gene/" ) # delete the output directory
n_significant_digits = 3

arr_data = ( df_expression * 10 ** n_significant_digits ).astype( int ).reset_index( drop = False ).values
# create the output folders
for path_folder in [ 
    f"{path_folder_wd}processed_data/web_application/base64/",
    f"{path_folder_wd}processed_data/web_application/gzipped/",
    f"{path_folder_wd}processed_data/expression/gene/",
] :
    os.makedirs( path_folder, exist_ok = True )
    
for arr in arr_data : # iterate each gene
    str_gene_symbol = arr[ 0 ]
    path_prefix = f"{path_folder_wd}processed_data/expression/gene/expr.multiplied_by_1{'0' * n_significant_digits}.{str_gene_symbol}"
    with gzip.open( f"{path_prefix}.tsv.gz", 'wb' ) as newfile :
        newfile.write( ( '\t'.join( [ str_gene_symbol ]+ list( map( str, arr[ 1 : ] ) ) ) + '\n' ).encode( ) ) # write an expression data of a file of the current gene_symbol
    Base64_Encode( f"{path_prefix}.tsv.gz", f"{path_prefix}.tsv.gz.base64.txt", header = ' ' ) # convert binary file into text using base64 encoding
    
''' assign each gene to different file if the size of file exceeds the given 'int_byte_file_size_limit' '''
# retrieve file size of base64 encoded expression data of each gene
df_file_base64 = GLOB_Retrive_Strings_in_Wildcards( f"{path_folder_wd}processed_data/expression/gene/expr.multiplied_by_1*.*.tsv.gz.base64.txt", retrieve_file_size = True )
df_file_base64.sort_values( 'wildcard_1', inplace = True ) # sort by gene_name

# assigne each gene to different file according to the given concatenated file size limit
int_byte_concatenated_file_size_limit = int( 90 * 2 ** 20 )
int_index_file, int_byte_accumulated = 0, 0
l_index_file = [ ]
for int_byte_file_size in df_file_base64.size_in_bytes.values :
    if int_byte_accumulated > int_byte_concatenated_file_size_limit :
        int_index_file += 1
        int_byte_accumulated = 0
    int_byte_accumulated += int_byte_file_size
    l_index_file.append( int_index_file )
df_file_base64[ 'index_file' ] = l_index_file

''' concatanate expression data for all genes and index the dataset '''
l_df_gene_index_byte = [ ]
str_significant_digits = df_file_base64.wildcard_0.values[ 0 ]
path_prefix_concatenated_file = f"{path_folder_wd}processed_data/web_application/base64/expr.multiplied_by_1{str_significant_digits}.tsv.gz.base64.concatenated.txt"
# for each separate concatenated file
for index_file in df_file_base64.index_file.unique( ) :
    # retrieve file information of base64 files belongining to the current concatenated file index
    df_file_base64_for_a_concatenated_file = PD_Select( df_file_base64, index_file = index_file )
    # concatanate base64 encoded files in the specified order
    OS_FILE_Combine_Files_in_order( df_file_base64_for_a_concatenated_file.path.values, f"{path_prefix_concatenated_file}.{index_file}", overwrite_existing_file = True )

    # write an index file describing the byte positions of each gene_symbol in the concatanated file
    int_byte_accumulated = 0
    l_l = [ ]
    for gene_symbol, size_in_bytes in df_file_base64_for_a_concatenated_file[ [ 'wildcard_1', 'size_in_bytes' ] ].values :
        l_l.append( [ gene_symbol, int_byte_accumulated, int_byte_accumulated + size_in_bytes ] ) # index_byte uses 0-based coordinates
        int_byte_accumulated += size_in_bytes # update accumulated number of bytes
    df_gene_index_byte = pd.DataFrame( l_l, columns = [ 'gene_symbol', 'index_byte_start', 'index_byte_end' ] )
    df_gene_index_byte[ 'index_file' ] = index_file 
    l_df_gene_index_byte.append( df_gene_index_byte )
# combine index for each index_file
name_prefix_index = f"expr.multiplied_by_1{str_significant_digits}.index"
df_gene_index_byte = pd.concat( l_df_gene_index_byte )
df_gene_index_byte.T.to_csv( f"{path_folder_wd}processed_data/web_application/gzipped/{name_prefix_index}.tsv.gz", sep = '\t', index = True, header = False )
Base64_Encode( f"{path_folder_wd}processed_data/web_application/gzipped/{name_prefix_index}.tsv.gz", f"{path_folder_wd}processed_data/web_application/base64/{name_prefix_index}.tsv.gz.base64.txt" ) # convert binary file into text using base64 encoding


# In[ ]:


# cells (identity, UMAP coordinates, and basic quality info) # 2021-07-14 20:44:16 
df_cells = df_cells.loc[ df_expression.columns.values ] # rearrange cell metatata in the order appears in the expression dataset
df_cells.index.name = 'id_cell'
df_cells.reset_index( inplace = True )
for col in [ 'UMAP_1', 'UMAP_2', 'tissue_UMAP_1', 'tissue_UMAP_2', 'percent_mt' ] : df_cells[ col ] = np.round( df_cells[ col ], 3 )
df_cells.T.to_csv( f"{path_folder_wd}processed_data/web_application/gzipped/cells.all.tsv.gz", sep = '\t', index = True, header = False )
Base64_Encode( f"{path_folder_wd}processed_data/web_application/gzipped/cells.all.tsv.gz", f"{path_folder_wd}processed_data/web_application/base64/cells.all.tsv.gz.base64.txt" ) # convert binary file into text using base64 encoding


# In[ ]:


# precompute Q1, Q2, Q3, min, max values of single cell RNA-Seq data for each tissue_type # 2020-11-28 19:28:35 
n_significant_digits = 3

df_expr = df_expression
arr_tissue = df_cells.tissue.values
l = list( )
for str_tissue in df_cells.tissue.unique( ) :
    df = df_expr.loc[ :, arr_tissue == str_tissue ]
    name_col = "{str_tissue} (n={n_cells})".format( str_tissue = str_tissue, n_cells = df.shape[ 1 ] )
    arr = np.vstack( [ df.min( axis = 1 ).values, df.quantile( q = 0.25, axis = 1, interpolation = 'nearest' ).values, df.median( axis = 1 ).values, df.quantile( q = 0.75, axis = 1, interpolation = 'nearest' ).values, df.max( axis = 1 ).values ] )
    df = pd.DataFrame( arr, columns = df.index.values, index = [ name_col + "|MIN", name_col + "|Q1", name_col + "|Q2", name_col + "|Q3", name_col + "|MAX" ] )
    l.append( df )
df_precomputed = pd.concat( l )
( df_precomputed * 10 ** n_significant_digits ).astype( int ).reset_index( drop = False ).to_csv( path_folder_wd + "processed_data/web_application/gzipped/all.box_precomputed_by_tissue.multiplied_by_1{n_significant_digits}.tsv.gz".format( n_significant_digits = "0" * n_significant_digits ), sep = '\t', index = False, header = False )


# In[ ]:


"""
Using expression data of each tissue data, perform correlation analysis
# 2024-02-28 17:14:49 
"""
# calculate correlation matrix (use a custom function utilizing precomputed metric for efficient calculation of correlation coefficient) # 2021-03-08 21:38:30 
def Correlation_Matrix_Pearsonr( arr, dtype = None ) :
    """ 
    # 2021-03-08 20:35:13 
    calculate correlation matrix of a given array (columns = cells/samples, rows = genes/observations/metrics)
    use precomputed mean-centered 'x' and 'y' divided by the matrix norm of itself and matrix multiplication to efficiently calculate correlation matrix of a gene expression matrix or other metric matrix
    """
    n = len( arr )
    if dtype is None : # default dtype np.float64
        dtype = np.float64
    # calculate mean-centered 'x' and 'y' by default
    # By using `astype(dtype)`, we ensure that the intermediate calculations
    # use at least 64 bit floating point.
    arr_mean = arr.mean( axis = 1, dtype = dtype )
    arr_mean_centered = arr - arr_mean.reshape( ( len( arr_mean ), 1 ) )
    
    # calculate matrix norm of mean-centered 'x' and 'y' by default
    # Unlike np.linalg.norm or the expression sqrt((xm*xm).sum()),
    # scipy.linalg.norm(xm) does not overflow if xm is, for example,
    # [-5e210, 5e210, 3e200, -3e200]
    # calculate mean-centered 'x' and 'y' divided by the matrix norm of itself by default
    
    arr_mean_centered_divided_by_norm = arr_mean_centered / np.vstack( list( scipy.linalg.norm( a ) for a in arr_mean_centered ) )
    arr_r = np.matmul( arr_mean_centered_divided_by_norm, arr_mean_centered_divided_by_norm.T )

    # Presumably, if abs(r) > 1, then it is only some small artifact of
    # floating point arithmetic.
    arr_r[ arr_r > 1 ] = 1
    arr_r[ arr_r < - 1 ] = - 1

    # As explained in the docstring, the p-value can be computed as
    #     p = 2*dist.cdf(-abs(r))
    # where dist is the beta distribution on [-1, 1] with shape parameters
    # a = b = n/2 - 1.  `special.btdtr` is the CDF for the beta distribution
    # on [0, 1].  To use it, we make the transformation  x = (r + 1)/2; the
    # shape parameters do not change.  Then -abs(r) used in `cdf(-abs(r))`
    # becomes x = (-abs(r) + 1)/2 = 0.5*(1 - abs(r)).  (r is cast to float64
    # to avoid a TypeError raised by btdtr when r is higher precision.)
    ab = n / 2 - 1
    arr_prob = np.array( list( scipy.special.btdtr( ab, ab, e ) for e in ( 1 - np.abs( arr_r.ravel( ) ) ) / 2 ), dtype = dtype ).reshape( arr_r.shape )

    return arr_r, arr_prob
def Pearsonr( x, y, xm = None, ym = None, normxm = None, normym = None, xm_divided_by_normxm = None, ym_divided_by_normym = None, dtype = None, return_intermediate_results = False ) :
    """ 
    assumes x and y are numpy arrays
    'xm' : precomputed mean-centered x values
    'ym' : precomputed mean-centered y values
    'normxm' : precomputed matrix norm of mean-centered x values
    'normym' : precomputed matrix norm of mean-centered x values
    'xm_divided_by_normxm' : mean-centered 'x' values divided by the matrix norm of itself
    'ym_divided_by_normym' : mean-centered 'y' values divided by the matrix norm of itself
    'return_intermediate_results' : return 'xm_divided_by_normxm' and 'ym_divided_by_normym' in addition to 'r' and 'prob'
    """
    n = len(x)
    
    if dtype is None : # default dtype is dtype of the array 'x'
        dtype = x.dtype
    # calculate mean-centered 'x' and 'y' by default
    # By using `astype(dtype)`, we ensure that the intermediate calculations
    # use at least 64 bit floating point.
    if xm is None and xm_divided_by_normxm is None :
        xmean = x.mean( dtype = dtype )
        xm = x.astype( dtype ) - xmean
    if ym is None and ym_divided_by_normym is None :
        ymean = y.mean( dtype = dtype )
        ym = y.astype( dtype ) - ymean
    # calculate matrix norm of mean-centered 'x' and 'y' by default
    # Unlike np.linalg.norm or the expression sqrt((xm*xm).sum()),
    # scipy.linalg.norm(xm) does not overflow if xm is, for example,
    # [-5e210, 5e210, 3e200, -3e200]
    if normxm is None and xm_divided_by_normxm is None :
        normxm = scipy.linalg.norm( xm )
    if normym is None and ym_divided_by_normym is None :
        normym = scipy.linalg.norm( ym )
    # calculate mean-centered 'x' and 'y' divided by the matrix norm of itself by default
    if xm_divided_by_normxm is None :
        xm_divided_by_normxm = xm / normxm
    if ym_divided_by_normym is None :
        ym_divided_by_normym = ym / normym

    r = np.dot( xm_divided_by_normxm, ym_divided_by_normym )

    # Presumably, if abs(r) > 1, then it is only some small artifact of
    # floating point arithmetic.
    r = max(min(r, 1.0), -1.0)

    # As explained in the docstring, the p-value can be computed as
    #     p = 2*dist.cdf(-abs(r))
    # where dist is the beta distribution on [-1, 1] with shape parameters
    # a = b = n/2 - 1.  `special.btdtr` is the CDF for the beta distribution
    # on [0, 1].  To use it, we make the transformation  x = (r + 1)/2; the
    # shape parameters do not change.  Then -abs(r) used in `cdf(-abs(r))`
    # becomes x = (-abs(r) + 1)/2 = 0.5*(1 - abs(r)).  (r is cast to float64
    # to avoid a TypeError raised by btdtr when r is higher precision.)
    ab = n / 2 - 1
    prob = 2 * scipy.special.btdtr( ab, ab, 0.5 * (1 - abs(np.float64( r ) ) ) )

    if return_intermediate_results :
        return r, prob, xm_divided_by_normxm, ym_divided_by_normym
    else :
        return r, prob
def perform_correl_analysis_and_save_results( df_expr, path_prefix_output : str ) :
    """ # 2023-06-28 19:01:06 
    perform correlation analysis and export outputs 
    """
    arr_r, arr_prob = Correlation_Matrix_Pearsonr( df_expr.values, dtype = np.float64 )
    arr_gene = df_expr.index.values
    df_r = pd.DataFrame( arr_r, index = arr_gene, columns = arr_gene )
    df_prob = pd.DataFrame( arr_prob, index = arr_gene, columns = arr_gene )
    df_r.to_csv( f"{path_prefix_output}.correlation_r.tsv.gz", sep = '\t' )
    df_prob.to_csv( f"{path_prefix_output}.correlation_prob.tsv.gz", sep = '\t' )  
    
path_folder_correl_output = f"{path_folder_wd}pipeline/correlation/"
os.makedirs( path_folder_correl_output, exist_ok = True ) # create the output folder

df_expr = df_expression
df_expr.index.name = '' # set default name for the index

""" for individual tissues """
arr_tissue = df_cells.tissue.values
for str_tissue in df_cells.tissue.unique( ) : # for each tissue
    perform_correl_analysis_and_save_results( df_expr.loc[ :, arr_tissue == str_tissue ], path_prefix_output = f"{path_folder_correl_output}tissue.{str_tissue}" )
    
""" for each combined dataset """
for name_combined, set_name_tissue in zip( [ 'mTissues6', 'Aging_and_Young_Kidney_Glom_Enrich' ],
[
    { 'Heart', 'Intestine', 'Kidney', 'Lung', 'Liver', 'Thymus' },
    { 'Glom-enriched (young mice)', 'Glom-enriched (old mice)', 'Tubule-enriched (young mice)', 'Tubule-enriched (old mice)', },
] ) :
    perform_correl_analysis_and_save_results( df_expr.loc[ :, list( e in set_name_tissue for e in arr_tissue ) ], path_prefix_output = f"{path_folder_correl_output}combined.{name_combined}" )
    
''' rename files so that the file name is compatible with URLs '''
# for e1, e2, e3, path_file in bk.GLOB_Retrive_Strings_in_Wildcards( f'{path_folder_wd}pipeline/correlation/tissue.* (*).correlation_*.tsv.gz' ).values :
#     os.rename( path_file, f"{path_folder_wd}pipeline/correlation/tissue.Kidney_{e1}_{e2}.correlation_{e3}.tsv.gz" )

''' simply renaming the correlation results, since the number of removed cells is 1% of all cells, and it will a lot of time '''
for o, n in zip( 
    [ 'Kidney_Glom_Old', 'Kidney_Tubule_Old', 'Kidney_Glom_Young', 'Kidney_Tubule_Young', ],
    [ 'Glom_enriched_young_mice', 'Glom_enriched_old_mice', 'Tubule_enriched_young_mice', 'Tubule_enriched_old_mice', ],
) :
    os.rename( f'{path_folder_wd}processed_data/web_application/base64/expr.correlation.tissue.{o}.index.tsv.gz.base64.txt', f'{path_folder_wd}processed_data/web_application/base64/expr.correlation.tissue.{n}.index.tsv.gz.base64.txt' )
    os.rename( f'{path_folder_wd}processed_data/web_application/base64/expr.correlation.tissue.{o}.tsv.gz.base64.concatanated.txt', f'{path_folder_wd}processed_data/web_application/base64/expr.correlation.tissue.{n}.tsv.gz.base64.concatanated.txt' )


# In[ ]:


''' encode correlation analysis results ''' # 2021-07-19 21:51:45 
# read correlation analysis result
os.makedirs( f"{path_folder_wd}pipeline/correlation/", exist_ok = True ) # create the output folder

def encode_correl_res( df_r, df_p, path_prefix_output : str, float_min_abs_r : float = 0.1, path_folder_temp : str = '/tmp/' ) :
    """ # 2023-06-30 00:04:23 
    encode correlation results
    
    float_min_abs_r = 0.1, # threshold minimum value of absolute value of correlation coefficients
    """
    # modified from WEB.Index_and_Base64_Encode
    # retrieve path_folder_wd
    path_folder_temp = os.path.abspath( path_folder_temp )
    if path_folder_temp[ -1 ] != '/' :
        path_folder_temp += '/'
    str_uuid = UUID( )
    path_folder_temp = f"{path_folder_temp}{str_uuid}/"
    os.makedirs( path_folder_temp, exist_ok = True )

    # write significant correlation results for each feature
    n_significant_digits_scientific_notation = 1
    n_significant_digits_typical_notation = 2

    arr_feature = df_r.index.values
    arr_r_all = df_r.values
    arr_p_all = df_p.values
    str_format_scientific_notation = "{:." + str( int( n_significant_digits_scientific_notation ) ) + "e}"
    str_format_typical_notation = '{:.' + str( int( n_significant_digits_typical_notation ) ) + 'f}'
    for index_feature in range( len( df_r ) ) :
        arr_r, arr_p = arr_r_all[ index_feature ], arr_p_all[ index_feature ]
        arr_mask = np.abs( arr_r ) > float_min_abs_r 
        arr_mask[ index_feature ] = False # ignore self-correlation (correlation coefficient is always 1)
        # since pandas operation is a bit slow, directoy compose tsv file and write as a gzip file
        str_content = ''
        str_content += "correlation_coefficient\t" + "\t".join( list( str_format_typical_notation.format( e ) for e in arr_r[ arr_mask ] ) ) + '\n'
        str_content += "probability\t" + "\t".join( list( str_format_scientific_notation.format( e ) for e in arr_p[ arr_mask ] ) ) + '\n'
        str_content += "gene_symbol\t" + "\t".join( arr_feature[ arr_mask ] ) + '\n'

        path_prefix_file = f"{path_folder_temp}{index_feature}"
        with gzip.open( f"{path_prefix_file}.tsv.gz", 'wb' ) as newfile :
            newfile.write( str_content.encode( ) )
        Base64_Encode( f"{path_prefix_file}.tsv.gz", f"{path_prefix_file}.tsv.gz.base64.txt", header = ' ' )

    # retrieve file size of base64 encoded chucks
    df_file_base64 = GLOB_Retrive_Strings_in_Wildcards( f"{path_folder_temp}*.tsv.gz.base64.txt", retrieve_file_size = True )
    df_file_base64.wildcard_0 = df_file_base64.wildcard_0.astype( int ) # retrieve integer index
    df_file_base64.sort_values( 'wildcard_0', inplace = True ) # sort by gene_name

    # concatanate base64 encoded files in the specified order
    OS_FILE_Combine_Files_in_order( df_file_base64.path.values, f"{path_prefix_output}.tsv.gz.base64.concatanated.txt", overwrite_existing_file = True )

    # write an index file describing the byte positions of each gene_symbol in the concatanated file
    int_byte_accumulated = 0
    l_l = [ ]
    for int_index, size_in_bytes in df_file_base64[ [ 'wildcard_0', 'size_in_bytes' ] ].values :
        l_l.append( [ arr_feature[ int_index ] ] + [ int_byte_accumulated, int_byte_accumulated + size_in_bytes ] ) # index_byte uses 0-based coordinates
        int_byte_accumulated += size_in_bytes # update accumulated number of bytes
    df_index_byte = pd.DataFrame( l_l, columns = [ 'gene_symbol' ] + [ 'index_byte_start', 'index_byte_end' ] )

    df_index_byte.T.to_csv( f"{path_folder_temp}index.tsv.gz", sep = '\t', index = True, header = False )
    Base64_Encode( f"{path_folder_temp}index.tsv.gz", f"{path_prefix_output}.index.tsv.gz.base64.txt" ) # convert binary file into text using base64 encoding

    # remove temporary folder
    shutil.rmtree( path_folder_temp )

for name_res in bk.GLOB_Retrive_Strings_in_Wildcards( f"{path_folder_wd}pipeline/correlation/*.correlation_r.tsv.gz" ).wildcard_0.values : # for each correlation result
    # read output results
    df_r = pd.read_csv( f"{path_folder_wd}pipeline/correlation/{name_res}.correlation_r.tsv.gz", sep = '\t', index_col = 0 )
    df_p = pd.read_csv( f"{path_folder_wd}pipeline/correlation/{name_res}.correlation_prob.tsv.gz", sep = '\t', index_col = 0 )
    print( f'reading data for {name_res} completed.' )
    
    # setting
    encode_correl_res( df_r, df_p, path_prefix_output = f"{path_folder_wd}processed_data/web_application/base64/expr.correlation.{name_res}", path_folder_temp = '/tmp/', float_min_abs_r = 0.2 )
    print( f'exporting output completed for {name_res} completed.' )


# In[ ]:


# encode gzipped binary files into base64-encoded gzipped files to enable transferring of the file as a plain text file.
# 2020-11-28 20:35:07 
for path_file_binary in glob.glob( path_folder_wd + "processed_data/web_application/gzipped/*.tsv.gz" ) :
    print( "processing {path_file_binary}".format( path_file_binary = path_file_binary ) )
    l = path_file_binary.rsplit( '/', 2 )
    l[ 1 ] = 'base64'
    path_file_binary_base64 = '/'.join( l ) + '.base64.txt'
    Base64_Encode( path_file_binary, path_file_binary_base64 ) # convert binary file into text using base64 encoding

