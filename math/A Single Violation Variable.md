### **Intro**

The computational constraints for FMP is too heavy, and we consider simplifying it by chaing the formulation for the constraint violation vector from one huge vector into just a single scalar. The new formulation is listed below: 

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
    \mathbf 1^T \lambda^{(r)} &\ge - \mathbf 1
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

---
### **LOG**


