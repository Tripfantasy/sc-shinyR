# sc-shinyR

**Note: Work in progress, features added and polished over time** - DM

A user-friendly single cell analysis app.

*Current test data is pbmc single cell assay via 10x genomics.*

Note: Make sure to update working directory in apprun.R before testing. Current upload limit is 10 gb. If uploading input files results in limit error, adjust options() in app.R to reflect the size of upload. eg., (10000*1024^2 == 10gb) 
Current state of app is WIP , expect errors and partial functionality as backend code is added. Functions are being called from a module.R , any functions should be added here to keep app.R clean.

### Functionality 
- Employs seurat functions for preprocessing and clustering, DEG analysis and visualization.
- Employs scType for cluster cell type annotation
- Custom thresholds for filtering, principal components, methods
- Opt-in to visualization inclusion in final report

### Roadmap 
- Basic app skeleton made [x]
- Preprocessing tab functionality [~]
- Annotation tab functionality []
- DEG tab functionality [] 
- Default values for configuration [] 
- KEGG / GO tab inclusion (necessary?) [] 
- Memory monitoring method [] 
- Variable input format functions []
- Raw data input to construct .Rds function []
- Aesthetic polishing []
- Interactive graphs in GUI []
- Graph 'snapshot' function to export to .pdf []
- Module optimization [] 
- Vignette example/walkthrough []
- Hover-over explanations of certain parameters. (E.g., Number of SCT features equates to the number of desired highly differentiated genes to be assessed in the SCT assay) Do we want these or just a good vignette/documentation. []
- Error messages (NaN/Inf/0 input)

#### Preprocessing (bug control) 
- redundancy in preprocessing function reloads object when it should already be loaded from load_obj
- hard-coded x as "nFeature_RNA" instead of dynamic x as pase0("nFeature_",assay) as temperary workaround to a "not found" error.
- hard coded graph.name argument to "RNA_snn" in FindClusters() as temporary workaround to "Graph.name not found in obj". This should be dynamic and based on output of FindNeighbors, nomenclature varied by default assay. 
