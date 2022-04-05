### **Master Problems Decision Variables and Constraints**

**Initial Start**: 

The master problems controls $w, \gamma$, $w$ is like: 

$$
w = \begin{bmatrix}
    x
    \\
    y
    \\ 
    z
\end{bmatrix}\quad 
x, y, z \in \mathbb \{0, 1\}^{|\mathcal G|\times |\mathcal T|}
$$

Therefore, $w$ should be able to be indexed by 3 sets $[s\in \{\text{x}, \text{y}, \text{z}\}, g \in \mathcal G, t\in \mathcal T]$. $\gamma$ is the objective, and it's just $\gamma \le M$, where $M$ is a vector ofthe same size, it's a huge positive vector to ensure that it's at least bounded. 

The parameter $\gamma$ **is a scaler** and it's directly put in as a scalar variable for the initial MP problem. 


**Adding Constraints as Feasibility Cut**



---
### **FMP: Limited FCP**

Parameters: $u, q$, continuous and discrete variables for secondary problem. All discrete variable are binary. Because the variables inside are not "uniform" (different index set), the variable will be implemented by a Dictionary mapping from the name of the variable to a vector representing that variable. 

Subvariables of $u$ are: 

$$
\begin{aligned}
    & p^{(t)}_n, c^{(t)}_n, (sr)_n^{(t)}, (regu)_n^{(t)}; (nsp)^{(t)}_n, (regd)_n^{(t)}, 
    (fft)_b^{(t)}, (frr)_b^{(b)}, (rfrr)_b^{(t)}, (tffr)_b^{(t)}, d^{(t)}_b
    \\
    & b \in \mathcal B; n \in \mathcal G; t \in \mathcal T
    \\
    & (inx)_{sg}^{(t)}
    \\
    & sg \in \mathcal N
    \\
    &\xi^+, \xi^- \in \mathbb R^{|\mathcal B|}, 
    \\
    & 
\end{aligned}
$$


Subvariables of $q$ are going to be:
$$
\begin{aligned}
    & \delta_{sg}^{(t)}, (x')_n^{(t)}, (y')_n^{(t)}, (z')_n^{(t)}
    \\
    & n\in \mathcal G', t\in \mathcal T, sg\in \mathcal N
    \\
    & \rho^+,\rho^-  \in \{0, 1\}^{|\mathcal B|}
\end{aligned}
$$

Both composite variables $u, q$ are going to be implemented via a dictionary mapping from the name of the variable as symbols for the key and maps to the instance for the variable in the code. 

**Constraints Matrix**


