"""
Simplified global Waddington-OT (gWOT) model.

This module implements the global optimization approach for trajectory inference
using entropic optimal transport between temporal marginals and data fitting terms.

Key components:
- PathsLoss: Regularization term based on OT between consecutive timepoints
- FitLoss: Data fitting term for a single timepoint
- TrajLoss: Combined loss function for the full trajectory
- optimize_model: Optimization routine using Adam

Students should implement the forward methods for each loss class.
"""

import numpy as np
import torch
from torch import nn
from torch.optim import Adam
from geomloss import SamplesLoss



# %%


# =========================
# 1. Path regularization
# =========================

class PathsLoss(nn.Module):
    """
    Compute the regularization term based on entropic optimal transport
    between consecutive timepoints.
    
    The path regularization enforces smoothness of the trajectory by penalizing
    large transport costs between consecutive marginals.
    
    Formula:
        L_reg = sum_{i=1}^{T-1} eps * Sinkhorn(P_{t_i}, P_{t_{i+1}})
    
    where eps = sigma * dt is the entropic regularization parameter.
    """

    def __init__(self, x_list, eps, device="cpu"):
        """
        Initialize path regularization loss.
        
        Parameters
        ----------
        x_list : list of torch.Tensor
            List of particle positions at each timepoint.
            x_list[i] has shape (n_particles, dim)
        eps : float
            Entropic regularization parameter (sigma * dt)
        device : str, optional
            Computation device ('cpu' or 'cuda')
        """
        super().__init__()
        self.x_list = x_list
        self.eps = torch.tensor(eps)
        self.device = device
        self.T = len(x_list)

        # TODO: Initialize the Sinkhorn loss from geomloss
        # Hint: Use SamplesLoss with loss="sinkhorn", p=2
        # For p=2, blur should be approximately sqrt(eps)
        # Use debias=False and scaling=0.95
        self.sinkhorn = SamplesLoss(loss = 'sinkhorn',p=2, blur = torch.sqrt(self.eps),scaling = 0.95,\
                                    debias=False)  
        
        

    def forward(self):
        """
        Compute the path regularization loss.
        
        Returns
        -------
        loss : torch.Tensor
            Total regularization loss (scalar)
        
        Implementation hints:
        - Loop over consecutive timepoints (t, t+1)
        - Compute Sinkhorn distance between x_list[t] and x_list[t+1]
        - Sum all distances
        - Multiply the total by self.eps
        """
        # TODO: Implement the forward pass
        total = 0.0
        # YOUR CODE HERE
        for i in range(self.T-1):
            total += self.sinkhorn(self.x_list[i],self.x_list[i+1])
        
        return total * self.eps


# =========================
# 2. Data-fitting loss 
# =========================

class FitLoss(nn.Module):
    """
    Compute the data fitting term for a single timepoint using a Gaussian
    likelihood formulation with log-sum-exp for numerical stability.
    
    The fitting loss measures how well the model particles explain the observed
    data by computing a negative log-likelihood under a Gaussian kernel.
    
    Formula:
        C_ij = -||x_i - y_j||^2 / (2 * sigma^2) + log(w_i)
        L_fit = -mean_j( log( sum_i exp(C_ij) ) )
    
    where:
    - x_i are model particles with weights w_i
    - y_j are observed particles
    - sigma is the standard deviation for the Gaussian kernel
    """

    def __init__(self, x, y, sigma, device="cpu"):
        """
        Initialize data fitting loss for one timepoint.
        
        Parameters
        ----------
        x : torch.Tensor, shape (n_model, dim)
            Model particle positions
        y : torch.Tensor, shape (n_obs, dim)
            Observed particle positions
        sigma : float
            Standard deviation for Gaussian likelihood
        device : str, optional
            Computation device
        """
        super().__init__()
        self.x = x
        self.y = y
        self.sigma = sigma
        self.device = device
        
    def forward(self):
        n_obs = self.y.shape[0]
        m = self.x.shape[0]
        

        # pairwise squared distances (m, n_obs)
        dist2 = ((self.x.unsqueeze(1) - self.y.unsqueeze(0))**2).sum(dim=2)
        
        # C_ij
        C = -dist2 / (2 * self.sigma**2)
        
        # second term
        lse = C - torch.logsumexp(C, dim=1).unsqueeze(1)
        term2 = - torch.logsumexp(lse, dim = 0).mean()
        
        # first term
        term1 = - torch.log(
            torch.tensor(n_obs / m, dtype=self.x.dtype, device=self.x.device)
            )
       
        return term1 + term2
        '''
    def forward(self):
        # pairwise squared distances (n_model, n_obs)
        dist2 = torch.cdist(self.x, self.y, p=2) ** 2
        
        # Gaussian log-kernel
        C = -dist2 / (2 * self.sigma**2)
        
        # add log weights (uniform weights)
        m = self.x.shape[0]
        log_w = -torch.log(torch.tensor(m, dtype=self.x.dtype, device=self.x.device))
        C = C + log_w   # broadcasts over columns
        
        # log-sum-exp over model particles (dim=0 → sum over i)
        log_likelihood = torch.logsumexp(C, dim=0)
        
        # mean over observations and negate
        return -log_likelihood.mean()
        '''

#FitLoss(torch.linspace(0,1,10),torch.linspace(0,1,20),0.1)

# =========================
# 3. Global TrajLoss (gWOT)
# =========================

class TrajLoss(nn.Module):
    """
    Combined loss function for global Waddington-OT (gWOT).
    
    This class combines the path regularization term and data fitting terms
    for all timepoints into a single objective function.
    
    Formula:
        L_total = lam_reg * L_reg + lam_fit * sum_{t=1}^T L_fit^(t)
    
    where:
    - L_reg is the path regularization (PathsLoss)
    - L_fit^(t) is the data fitting term at time t (FitLoss)
    - lam_reg and lam_fit are weighting coefficients
    """

    def __init__(
        self,
        x0_list,
        obs_list,
        lam_reg,
        lam_fit,
        eps,
        sigma_fit=1.0,
        device="cpu",
    ):
        """
        Initialize global trajectory loss (gWOT).
        
        Parameters
        ----------
        x0_list : list of torch.Tensor
            Initial particle positions for each timepoint.
            Each tensor has shape (n_particles, dim)
        obs_list : list of torch.Tensor
            Observed particles at each timepoint.
            Each tensor has shape (n_obs_t, dim)
        lam_reg : float
            Regularization weight (controls smoothness)
        lam_fit : float
            Data fitting weight (controls fidelity to observations)
        eps : float
            Entropic regularization parameter (sigma * dt)
        sigma_fit : float, optional
            Standard deviation for data fitting Gaussian kernel
        device : str, optional
            Computation device
        
        Notes
        -----
        The particle positions (x0_list) should be stored as nn.Parameter
        so that PyTorch can optimize them during training.
        """
        super().__init__()

        self.device = device
        self.T = len(x0_list)

        # TODO: Convert x0_list to nn.ParameterList
        # Hint: Use nn.ParameterList with nn.Parameter for each tensor
        # Make sure tensors are on the correct device and dtype (float32)
        self.x = nn.ParameterList([
            nn.Parameter(x.to(device=self.device, dtype=torch.float32))
            for x in x0_list
        ])
        # Store observations
        self.obs_list = [
            obs.to(device=self.device, dtype=torch.float32).contiguous()
            for obs in obs_list
        ]

        # Store hyperparameters
        self.lam_reg = lam_reg
        self.lam_fit = lam_fit
        self.eps = eps
        self.sigma_fit = sigma_fit

    def forward(self):
        """
        Compute total loss.
        
        Returns
        -------
        loss : torch.Tensor
            Total loss (scalar)
        
        Implementation hints:
        - Create PathsLoss with self.x and self.eps, call it to get reg_loss
        - Loop over timepoints to compute fit_loss:
            * Create FitLoss with model particles, observations, sigma_fit
            * Accumulate the fit losses
        - You may want to scale fit_loss by the number of observations for balance
        - Combine: lam_reg * reg_loss + lam_fit * fit_loss
        """
        # TODO: Implement the forward pass
        # YOUR CODE HERE
        
        reg_loss = PathsLoss(self.x, self.eps)()
        fit_loss = torch.tensor(0.0)
        for i in range(self.T):
            #print(FitLoss(x = self.x[i], y = self.obs_list[i], sigma = self.sigma_fit))
            fit_loss += FitLoss(x = self.x[i], y = self.obs_list[i], sigma = self.sigma_fit)()
        #fit_loss /= self.T

        return self.lam_reg * reg_loss + self.lam_fit * fit_loss
    
# =========================
# 4. Optimization with Adam
# =========================

def optimize_model(
    traj_model,
    n_epochs=1000,
    lr=1e-2,
    print_every=100,
):
    """
    Optimize the trajectory model using Adam optimizer.
    
    This function performs gradient descent on the particle positions to minimize
    the gWOT objective function.
    
    Parameters
    ----------
    traj_model : TrajLoss
        The trajectory loss model to optimize
    n_epochs : int, optional
        Number of optimization iterations
    lr : float, optional
        Learning rate for Adam optimizer
    print_every : int, optional
        Print loss every N iterations
    
    Returns
    -------
    history : list of float
        Loss values at each iteration
    best_positions : list of np.ndarray
        Optimized particle positions for each timepoint
    
    Notes
    -----
    This function uses the Adam optimizer which adapts the learning rate
    for each parameter. It typically converges faster than vanilla SGD.
    """
    
    # Initialize Adam optimizer on model parameters
    optimizer = Adam(traj_model.parameters(), lr=lr)
    history = []

    # Initial gradient computation (warm-up)
    traj_model.zero_grad(set_to_none=True)
    loss = traj_model()
    loss.backward()

    # Main optimization loop
    for epoch in range(n_epochs):
        # Zero gradients from previous iteration
        optimizer.zero_grad()
        
        # Forward pass: compute loss
        loss = traj_model()
        
        # Backward pass: compute gradients
        loss.backward()
        
        # Update parameters
        optimizer.step()

        # Store loss and best positions
        with torch.no_grad():
            loss_val = loss.item()
            history.append(loss_val)
            
            # Extract current particle positions
            best_positions = [
                x_t.detach().cpu().numpy().copy() for x_t in traj_model.x
            ]

        # Print progress
        if print_every > 0 and (epoch + 1) % print_every == 0:
            print(f"[Adam] Epoch {epoch+1}/{n_epochs}, loss = {loss_val:.4f}")

    return history, best_positions
