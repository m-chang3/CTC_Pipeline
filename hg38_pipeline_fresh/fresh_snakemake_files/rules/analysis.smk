rule create_Seurat:
     input:
        "data/gene_metadata.tsv",
        "data/counts/raw_counts_.filt.tsv"
     output:
       "plots/seurat/QC_vlnplot0.pdf",
       "plots/seurat/QC_scatter_percentMT0.pdf",
       "plots/seurat/QC_scatter_nfeature0.pdf",
       "plots/seurat/Heatmap10.pdf",
       "plots/seurat/Elbowplot0.pdf",
       "plots/seurat/UMAP0.pdf",
       "plots/seurat/top10_vlnplot0.pdf",
       "data/seurat/SeuratObject.rds"          

     params:
        coverage_threshold = config["coverage_threshold"],
        features_threshold = config["features_threshold"],
        top50_threshold = config["top50_threshold"],
        Count_upperQuantile = config["Count_upperQuantile"],
        Feature_lowerQuantile = config["Feature_lowerQuantile"],
        Feature_upperQuantile = config["Feature_upperQuantile"],
        percentMT_upperQuantile = config["percentMT_upperQuantile"],
	integrateTF = config["integrateTF"],
	clusterRes = config["clusterRes"]
     conda:
        "../envs/updated_seurat.yaml"        	
     shell:
        "Rscript scripts/create_seurat_matt.R --covThresh={params.coverage_threshold} --featThresh={params.features_threshold} --top50={params.top50_threshold} --featureLQ={params.Feature_lowerQuantile} --featureUQ={params.Feature_upperQuantile} --CountUQ={params.Count_upperQuantile} --MTUQ={params.percentMT_upperQuantile} --integrate={params.integrateTF} --res={params.clusterRes}"
