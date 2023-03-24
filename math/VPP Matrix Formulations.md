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


### **Matrix $B$, Vector $\bar{w}$ | First Stage Binary**
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


### **Matrix $C$, Vector $u$ | Second Stage Continuous**

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
        (dr)_t
        \\
        \mu_{l}^{(t)}
    \end{bmatrix}
\end{aligned}
$$


### **Matrix $G$, Vector $q$ | Second Stage Binary**

Secondary discrete variables, they are for the secondary constraints and the secondary generators, they are packed into the $q$ vector. 

$$
\begin{aligned}
    q = 
    \begin{bmatrix}
        (\sigma^+_s)^{(t)}
        \\
        (\sigma^-_s)^{(t)}
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


--- 
### **First Stage Constraints**

$$
\begin{aligned}
    & x_{n}^{t} + z_{n}^{t} \leq 1, \forall n \in \mathcal{G}, \forall t,
    \\
    & y_{n}^{t} - y_{n}^{t-1} = x_{n}^{t} - z_{n}^{t},
    \forall n \in \mathcal{G}, \forall t, 
    \\
    & {\sum_{\tau = t - T_{n}^{\text{T},\text{minu}} + 1}^{t} x_{n}^{\tau} \leq y_{n}^{t},}
    \quad 
    { \forall n \in \mathcal{G}, \forall t \in \{T_{n}^{\text{ T},\text{ minu}}, \cdots, T_{end}\}}, \\
    &{\sum_{\tau = t - T_{n}^{\text{T},\text{mind}} + 1}^{t} z_{n}^{\tau} \leq 1 - y_{n}^{t},} 
    \quad 
    {\forall n \in \mathcal{G}, \forall t \in \{T_{n}^{\text{T},\text{ mind}}, \cdots,T_{end}\}}
\end{aligned}
$$


### **Fuel Constraints**


$$
\begin{aligned}
    & 
    \underbrace{
        \sum_{t \in \mathcal T} 
        \sum_{n \in {\mathcal{G}}}
        \text{SU}_{n} x_{n}^{(t)} 
        + 
        \sum_{t \in \mathcal T} 
        \sum_{n \in {\mathcal{G}}}\text{SD}_{n} z_{n}^{(t)} 
    }_{\text{for }B}
    + 
    \underbrace{
        \sum_{t \in \mathcal T} 
        \sum_{n \in {\mathcal{G}}}c_{n}^{t}
        +
        \sum_{t \in \mathcal T} 
        \sum_{l=1}^L\rho_{l} \mu_{l, t}
    }_{\text{for }   C}
    \leq  \phi
    \\
    & 
    \underbrace{\alpha_n^{(k)}p_n^{(t)}- c_n^{(t)}}_{\text{for }C}
    + 
    \underbrace{\beta_n^{(k)}y_n^{(t)} }_{\text{for }B}
    \le 0
\end{aligned}
$$


### **Capacity constraints:** 
$$
\begin{aligned}
    &
    \underbrace{p_n^{(t)} - p^{(t - 1)}_n }_{\text{for }C}
    \underbrace{- \text{RU}_ny_n^{(t - 1)}
    - \overline{\text{RU}}_n x_n^{(t)} }_{\text{for }B}
    \le 0
    \\&
    \underbrace{p^{(t - 1)}_n - p^{(t)}_n }_{\text{for }C}
    \underbrace{- \text{RD}_ny_n^{(t)} -
    \overline{\text{RD}}_nz_n^{(t)} }_{\text{for }B}
    \le 0
    \\&
    \underbrace{p_n^{(t)}}_{\text{for }C}
    \underbrace{- P^{\max}_ny_n^{(t)}}_{\text{for } B} \le 0
    \\&
    \underbrace{- p_n^{(t)}}_{\text{for }C} 
    + \underbrace{y_n^{(t)}P_n^{\min}}_{\text{for } B} \le 0
\end{aligned}
$$


### **Demand Response Constraints**

$$
\begin{aligned}
    & 
    - \mu_{l, t} \le 0 &\quad \forall t \in [T], l \in [L]
    \\
    &
    \mu_{l,t}\leq R_{t,l}-R_{t,l-1} &\quad \forall t \in [T], l \in [L]
    \\
    &
    (dr)_{t} - \sum_{l=1}^L \mu_{l,t} = 0
    \\
    & 
    (dr)_t\leq DR_{max}
\end{aligned}
$$
constraints are all are in the same matrix: $C$. 


### **Battery Constraints**

$$
\begin{aligned}
    & 
    \underbrace{h_s^{1} - h_s^{0} = 0}_{\text{for C}}, \forall s \in S \quad 
    \\
    \text{for }C
    & 
    \begin{cases}
        h_s^{t} - h_s^{t - 1} - \nu_s^{+}(g_s^+)^{t-1} + \nu_s^{-}(g_s^-)^{t-1} = 0, 
        & 
        \hspace{3mm} 
        \forall t=2,\cdots,T, \forall s 
        \\
        h_s^t \leq \bar{H}_s & \hspace{3mm} \forall t 
    \end{cases}
    \\
    & 
    \underbrace{(g^+)_s^{(t)}}_{\text{for }C}
    -
    \underbrace{(\sigma_s^+)^{(t)}\bar{G}_s^{+}}_{\text{for }G}
    \leq 0  
    \hspace{3mm} \forall t,s
    \\
    & \underbrace{(g^-)_s^{(t)}}_{\text{for } C}
    -
    \underbrace{(\sigma_s^-)^{(t)}\bar{G}_s^{-}}_{\text{for }G}\le 0 \hspace{3mm} \forall t,s
    \\
    &
    \underbrace{
        (\sigma_s^+)^{(t)} +
        (\sigma_s^-)^{(t)}
    }_{\text{for } G}=1 \hspace{3mm} \forall t,s
\end{aligned}
$$

### **Demand Balance Constraints**: 

$$
\begin{aligned}
    \sum_{n \in \mathcal G}
    p_{n}^{t} 
    +
    \left(
    \sum_{s\in \mathcal{S}}
        (g^-)_s^{t}-(g^+)_s^{(t)}
    \right)
    +
    (dr)_t
    =
    d_{t}
    \quad \forall t
\end{aligned}
$$

They are all for the matrix $C$. 
