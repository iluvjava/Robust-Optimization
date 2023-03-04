### **Over All Problem**

Robust statement: 

$$
\max (\gamma) \; \text{s.t: }
\\
Aw \leq b
\\
Bw+Cu+Gq \leq Hd+h
\\
\forall d \in [\hat{d}-\gamma\bm{1}, \hat{d}+\gamma\bm{1}]
\\
\begin{cases}
    w, &\text{binary}
    \\
    u\ge \mathbf 0, &\text{continuous}
    \\
    q, &\text{binary}
\end{cases}
$$


### **Matrix $B$, Vector $\bar{w}$**
Note, the vector $\bar{w}$ is not a decision variable, it's a fixed value passed from the master problem determining the first stage decision variables. 

$$
\begin{aligned}
    \overline{w} = \begin{bmatrix}
        x_n^{(t)} 
        \\
        y_n^{(t)}
        \\
        z_n^{(t)}
    \end{bmatrix}
\end{aligned}
$$


### **Matrix $C$, Vector $u$**

Secondary continuous variables and their constraints. 

$$
u = 
\begin{aligned}
    \begin{bmatrix}
        c^{(t)}_n
        \\
        p^{(t)}_n
        \\
        h^{(t)}_s
        \\
        (g^+)^{(t)}_s
        \\
        (g^-)^{(t)}_s
        \\
        dr_t
        \\
        \mu_{l, t}
    \end{bmatrix}
\end{aligned}
$$


### **Matrix $G$, Vector $q$**

Secondary discrete variables, they are for the secondary constraints and the secondary generators, they are packed into the $q$ vector. 

$$
\begin{aligned}
    q = 
    \begin{bmatrix}
        \sigma^+_t
        \\
        \sigma^-_t
    \end{bmatrix}
\end{aligned}
$$


### **Matrix $H$, Vector $d$**

The demand vector for the feasibility problem: 

$$
\begin{aligned}
    \begin{bmatrix}
        d_t
    \end{bmatrix}
\end{aligned}
$$

Index by the time. 