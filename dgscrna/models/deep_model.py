"""
Deep learning model for DGscRNA package
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Optional

class DeepModel(nn.Module):
    """
    Deep neural network for cell type annotation refinement
    """
    
    def __init__(self, input_dim: int, num_classes: int, hidden_dims: Optional[list] = None):
        """
        Initialize the deep model
        
        Parameters
        ----------
        input_dim : int
            Input dimension (number of genes)
        num_classes : int
            Number of cell type classes
        hidden_dims : list, optional
            List of hidden layer dimensions. Default: [256, 128]
        """
        super().__init__()
        
        if hidden_dims is None:
            hidden_dims = [256, 128]
        
        # Build layers
        layers = []
        prev_dim = input_dim
        
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.LeakyReLU(),
                nn.Dropout(0.3)
            ])
            prev_dim = hidden_dim
        
        # Output layer
        layers.append(nn.Linear(prev_dim, num_classes))
        
        self.network = nn.Sequential(*layers)
        
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass
        
        Parameters
        ----------
        x : torch.Tensor
            Input tensor of shape (batch_size, input_dim)
            
        Returns
        -------
        torch.Tensor
            Output tensor of shape (batch_size, num_classes) - raw logits
        """
        return self.network(x)
    
    def predict(self, x: torch.Tensor) -> torch.Tensor:
        """
        Make predictions
        
        Parameters
        ----------
        x : torch.Tensor
            Input tensor
            
        Returns
        -------
        torch.Tensor
            Predicted class indices
        """
        with torch.no_grad():
            outputs = self.forward(x)
            return torch.argmax(outputs, dim=1)
    
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """
        Get prediction probabilities
        
        Parameters
        ----------
        x : torch.Tensor
            Input tensor
            
        Returns
        -------
        torch.Tensor
            Prediction probabilities
        """
        with torch.no_grad():
            logits = self.forward(x)
            return F.softmax(logits, dim=1) 