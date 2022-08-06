#运行程序

#文件必须是.txt结尾格式
cellphonedb method statistical_analysis --threads 20 --counts-data ensembl --output-path test_output meta.txt exprMatrix_human3.txt

##气泡图
cellphonedb plot dot_plot --pvalues-path test_output/pvalues.txt  --means-path test_output/means.txt --output-path test_output  --output-name test.dotplot.pdf

##热图
cellphonedb plot heatmap_plot --pvalues-path test_output/pvalues.txt --output-path test_output --pvalue 0.05 --count-name test.heatmap_count.pdf --log-name test.heatmap_log_count.pdf --count-network-name test.count_network.txt --interaction-count-name test.interaction_count.txt meta.txt
#
