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
    \lambda^{(r)} &\le \mathbf 0 \quad &\forall r = 1, \cdots, k
    \\
     \hat d - \vec{\gamma} &\le d \le \hat d + \vec \gamma \quad &\forall i = 1, \cdots, k
\end{aligned}
$$

Defines 

- FMPH1: $d$ is a constant and $\lambda$ is a decision variables. 
- FMPH2: $\lambda$ is a constant from FMPH1 and $d$ is a continuous decision variables. 

**Warm Start**: 

To start with an already feasible solution from the binlinear heuristic on the FMP problem, check out: [this](https://jump.dev/JuMP.jl/stable/tutorials/conic/start_values/#Primal-and-dual-warm-starts). 