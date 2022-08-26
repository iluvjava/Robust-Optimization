### **Intro**

The originally conceive the cut: 

$$
    Cu_*^{(k + 1)} + Gq_*^{(k + 1)} \le H\hat d + H\Gamma ((\rho_*^+)^{(k + 1)} - (\rho_*^-)^{(k + 1)}) + h - Bw
$$

where $u_*, q_*$ are the latest decision variables fixed from the later iterations of FSP, CCGA inner for loop. To avoid the violation vector $v$ causing problem here, we make $u^{(i)}, k^{(i)}$ a decision variable of its own, indexed by $i$ the outter CCGA iterations for the cut on the master problem. 

When FMP, FSP returns positive $v$, the above equation cannot be satisfied with the given $w$ for all $u, q$, hence giving the control to $u, q$ will eliminate the choice of $w\bar$ which initialized the CCGA inner forloop iterations. 
