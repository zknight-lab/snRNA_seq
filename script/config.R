###################
###Set Work Path###
###################


workpath='d:/snRNA_seq/'
datapath='d:/snRNA_seq/Data/'
RData='d:/snRNA_seq/RData/'
run_1='run_1/filtered_feature_bc_matrix/'
run_custom='run_custom/filtered_feature_bc_matrix/'
run_tdtomato='run_tdtomato/filtered_feature_bc_matrix/'
run_VDB_1='run_VDB_1/filtered_feature_bc_matrix/'
output='d:/snRNA_seq/output/'
cluster='cluster/'
script='script/'

setwd(workpath)

###init_config###----

normalization.method = 'LogNormalize'
selection.method = 'vst'
project <- 'forebrain'
dim = 30
seed.use= NULL
scale.factor = 10000
nfeatures = 2000


filename = 'Zhuang'