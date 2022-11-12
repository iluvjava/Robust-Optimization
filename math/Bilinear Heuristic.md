### **Bilinear Heuristic**
$$
\begin{aligned}
     \varphi(\bar w, \bar \gamma&) = \max_{
            \substack{
            u^{r}, v^{r}, q^{r}, \lambda^{(r)}\ge \mathbf 0\\
            \forall r = 1, \cdots k
            \\
            d
            }
        }\eta
    \\
    \eta &\le v^{(r)} \quad &\forall r = 1, \cdots, k
    \\
    Cu^{(r)} - \mathbf 1 v^{(r)} &\le 
    Hd + h - B \bar w - Gq^{(r)}
    \quad 
    &\forall r = 1,\cdots, k
    \\
    (\lambda^{(r)})^T C &\le \mathbf 0 \quad &\forall i = 1, \cdots, k
    \\
    \mathbf 1^T \lambda^{(r)} &\ge -1
    \quad 
    &\forall i = 1, \cdots, k
    \\
    (\lambda^{(r)})^T(Hd + h - B \bar w - Gq^{(r)}) 
    &= v^{(r)} 
    \quad 
    &\forall r = 1,\cdots, k
    \\
    \lambda^{(r)} &\le \mathbf 0 \quad &\forall r = 1, \cdots k
    \\
     \hat d - \vec{\gamma} &\le d \le \hat d + \vec \gamma \quad &\forall i = 1, \cdots, k
\end{aligned}
$$


1. Make the model, set $d$ to be the boundary of their interval, randomly. 
* while wihtin a certain iteration limit: `N`  and `eta == 0`
   1. Solve the system (keep an eye for $λ$). 
   2. If objective is positive then give the results to FSP. 
   3. Set $λ^{(r)}$ to be a constant the equal so to the solution from previous solve. Change $d$ back to a continuous decision variables. 
   4. Solve it again with $\lambda^{(r)}$ 
   5. Project $d$ to the boundary.  

**Warm Start**: 

To start with an already feasible solution from the binlinear heuristic on the FMP problem, check out: [this](https://jump.dev/JuMP.jl/stable/tutorials/conic/start_values/#Primal-and-dual-warm-starts). 