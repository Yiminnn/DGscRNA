## **DGscRNA Algorithm Pipeline Overview**

DGscRNA (Density-based Gene expression scoring for single-cell RNA-seq with Deep Learning) is a two-stage cell type annotation method that integrates marker-based density scoring with neural network classification.

---

## **Stage 1: Preprocessing & Clustering**

### **1.1 Data Preprocessing**
```r
# Quality control and normalization
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern="^MT-")
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt < 15)
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize")
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
```

### **1.2 Dimensionality Reduction & Clustering**
- **PCA**: Principal component analysis for initial dimensionality reduction
- **UMAP**: Non-linear dimensionality reduction for visualization
- **Multiple Clustering**: Supports various clustering algorithms:
  - Seurat (graph-based clustering)
  - HDBSCAN (density-based clustering)
  - Monocle3 (trajectory-based clustering)
  - scCCESS (consensus clustering)

---

## **Stage 2: Density-Based Scoring**

### **2.1 Marker Gene Loading**
The algorithm loads cell type marker genes from multiple databases:
- **CellMarker 2.0**: Curated cell type markers
- **Human Protein Atlas**: Tissue-specific markers
- **PanglaoDB**: Single-cell marker database
- **Custom marker sets**: User-defined markers

### **2.2 Density Score Calculation**
For each cluster and cell type combination:

```r
density_score <- function(integrated_data, ct_markers, DEG_markers_set, clustering, cutoffs) {
    # For each cluster:
    curr_markers <- DEG_markers[DEG_markers$cluster == cluster_i, ]
    curr_markers <- curr_markers[curr_markers$avg_log2FC > 1, ] # Upregulated genes
    
    # For each cell type marker set:
    intersecting_genes <- intersect(curr_cell_type_set, top_genes)
    
    # Calculate density score:
    score = sum(avg_log2FC_of_intersecting_genes) / length(cell_type_marker_set)
}
```

### **2.3 Multi-Threshold Scoring**
Three cutoff strategies are applied:
- **`0.5`**: Hard threshold (score ≥ 0.5)
- **`mean`**: Dynamic threshold (score ≥ mean of all scores)
- **`none`**: No threshold (accepts highest score)

### **2.4 Initial Annotation**
```r
# Assign cell types based on maximum density scores
if (max_score < threshold) {
    annotation <- "Undecided"
} else {
    annotation <- cell_type_with_max_score
}
```

---

## **Stage 3: Deep Learning Refinement**

### **3.1 Data Preparation**
```python
# Separate cells into training and prediction sets
agreed_cells = cells_with_confident_annotations  # Training data
disagreed_cells = cells_labeled_as_"Undecided"   # Prediction targets

X_train = gene_expression_matrix[agreed_cells, :]
y_train = cell_type_labels[agreed_cells]
X_test = gene_expression_matrix[disagreed_cells, :]
```

### **3.2 Neural Network Architecture**
```python
class DeepModel(nn.Module):
    def __init__(self, input_shape, class_shape):
        self.fc1 = nn.Linear(input_shape, 256)    # Input: 2000 genes
        self.fc2 = nn.Linear(256, 128)            # Hidden layer
        self.fc3 = nn.Linear(128, class_shape)    # Output: cell types
        self.activation = nn.LeakyReLU()
        self.output_activation = nn.Softmax()
```

### **3.3 Training Process**
```python
for epoch in range(EPOCHS):
    for batch in train_data:
        preds = model(samples)
        loss = criterion(preds, labels)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
```

### **3.4 Confidence-Based Prediction**
```python
# Apply trained model to "Undecided" cells
probabilities = model(X_test)
max_probs = np.max(probabilities, axis=1)

# High confidence predictions (≥90%) are accepted
confident_predictions = predictions[max_probs >= 0.90]
low_confidence = predictions[max_probs < 0.90]  # Remain "Unknown"
```

---

## **Stage 4: Integration & Output**

### **4.1 Final Annotation**
The algorithm produces multiple annotation columns:
- `{clustering}_{marker_source}_{celltype}_{cutoff}`: Density-based annotations
- `{clustering}_{marker_source}_{celltype}_{cutoff}_DGscRNA`: Deep learning refined annotations

### **4.2 Quality Metrics**
- **Accuracy**: Overall prediction accuracy
- **F1-Score**: Balanced precision and recall
- **Precision/Recall**: Per-class performance metrics
- **Confidence scores**: Prediction reliability

---

## **Key Algorithm Advantages**

1. **Multi-Source Integration**: Combines multiple marker databases
2. **Adaptive Thresholding**: Multiple cutoff strategies for robustness
3. **Deep Learning Refinement**: Improves uncertain predictions
4. **Confidence Scoring**: Provides prediction reliability estimates
5. **Scalable**: Works with various clustering methods and datasets

---

## **Workflow Summary**

```
Single-cell Data → Preprocessing → Clustering → Marker Loading
                                      ↓
Density Scoring → Multi-threshold Evaluation → Initial Annotations
                                      ↓
"Undecided" Cells → Deep Learning Training → Confidence-based Prediction
                                      ↓
Final Annotations with Confidence Scores
```

This pipeline effectively combines the interpretability of marker-based methods with the power of deep learning to provide robust, confident cell type annotations for single-cell RNA-seq data.