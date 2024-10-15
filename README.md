# sc-shinyR (WIP)

**Note: Project will be re-visited at some point in the future.** - DM 2024

A user-friendly single cell analysis app.

*Current test data is pbmc single cell assay via 10x genomics.*

Note: Make sure to update working directory in apprun.R before testing. Current upload limit is 10 gb. If uploading input files results in limit error, adjust options() in app.R to reflect the size of upload. eg., (10000*1024^2 == 10gb) 
Current state of app is WIP , expect errors and partial functionality as backend code is added. Functions are being called from a module.R , any functions should be added here to keep app.R clean.

### Functionality 
- Employs seurat functions for preprocessing and clustering, DEG analysis and visualization.
- Employs scType for cluster cell type annotation
- Custom thresholds for filtering, principal components, methods
- Opt-in to visualization inclusion in final report

