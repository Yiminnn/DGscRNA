"""
scanpy-dgscrna: A scanpy extension for automated cell type annotation using density scoring and deep learning.

Example usage:

```python
import scanpy as sc
import scanpy_dgscrna as dgsc

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Run the complete workflow
dgsc.run_dgscrna_workflow(
    adata, 
    marker_folder='path/to/markers',
    cluster_columns=['leiden', 'seurat_clusters']
)

# Or run individual steps
marker_sources = dgsc.load_markers('path/to/markers')
dgsc.density_score(adata, marker_sources, ['leiden'])
dgsc.dgscrna_annotate(adata, 'leiden_CellMarker_Bcell_0.5')
```
"""

from .tools import *
from .plotting import *
from .preprocessing import *
