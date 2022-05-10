### **Original Form of the Problem**
How is the original form of the problem and how it's phrased??? What are these matrices part of??? Here is what the matrices are trying to make: 

$$
\begin{aligned}
    & 
    \max_d\min_{q\in Q}\min_{u\ge \mathbf 0, v\mathbf \ge \mathbf 0}
    \langle \mathbf 1, v\rangle
    \\
    & 
    \hat{d} - \bar{\gamma} \le d \le \hat{d} - \bar\gamma
    \\
    & 
    Cu - v + B\bar{w} + Gq\le d + h
\end{aligned}
$$

Maximizes the demand $d$ to keep the problem feasible for all $q\in Q$ for all secondary continuous variable $u$, and slack variable $v\ge \mathbf 0$. 

---
### **Enumerating Sets**

---
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



----
### **Matrix $C$, Vector $u$**

Secondary continuous variables and their constraints. 

$$
u = 
\begin{aligned}
    \begin{bmatrix}
        c^{(t)}_n
        \\
        (c')^{(t)}_m
        \\
        p^{(t)}_n
        \\
        (p')^{(t)}_m
        \\
        (regu)^{(t)}_n
        \\
        (regu')^{(t)}_m
        \\
        (regd)^{(t)}_n
        \\
        (regd')^{(t)}_m
        \\
        (sr)^{(t)}_n
        \\
        (sr')^{(t)}_m
        \\
        h^{(t)}_s
        \\
        (g^+)^{(t)}_s
        \\
        (g^-)^{(t)}_s
        \\
        (nsp)^{(t)}_n
        \\
        (nsp')^{(t)}_m
    \end{bmatrix}
\end{aligned}
$$

Variable range: 

$$
\begin{aligned}
    p_n, regu, regd, sr \ge \mathbf 0
\end{aligned}
$$

---
### **Dangling Variables**
Decision variables that are not involved with any matrices but they still appear in the original master problem formulations.  

* $d$ is the demand variables. Directly affecting the RHS of all the constraints. Only appears on the RHS of the demand balance constraints set. 


---
### **Matrix $G$, Vector $q$**

Secondary discrete variables, they are for the secondary constraints and the secondary generators, they are packed into the $q$ vector. 

$$
\begin{aligned}
    q = 
    \begin{bmatrix}
        (x')^{(t)}_m 
        \\
        (y')^{(t)}_m
        \\
        (z')^{(t)}_m
    \end{bmatrix}
\end{aligned}
$$

---
### **Constraints Classifications**

Fuel Constraints: 

$$
\begin{aligned}
    & 
    \underbrace{\sum_{t \in \mathcal T}^{}\sum_{m\in \mathcal G'}^{}
    \left(
        \text{SU}_m(x')^{(t)}_m
    \right)
    + 
    \sum_{t \in \mathcal T}^{}\sum_{m\in \mathcal G'}^{}
    \left(
        \text{SD}_m(z')^{(t)}_m
    \right)}_{\text{for }G}
    \\
    &
    +
    \underbrace{\sum_{t \in \mathcal T}^{}\sum_{m\in \mathcal G'}^{}
        (c')^{(t)}_m
    +
    \sum_{t \in \mathcal T}^{}
    \sum_{n\in \mathcal G}^{}
    c^{(t)}_n }_{\text{for }C}
    \\&
    +
    \underbrace{\sum_{t \in \mathcal T}^{}
    \sum_{n\in \mathcal G}^{}
    \left(
        \text{SU}_nx_n^{(t)} + \text{SD}_nz_n^{(t)} \tag{6} 
    \right)}_{\text{for } B}
    \le
    \phi
\end{aligned}
$$


$$
\begin{aligned}
    & 
    \underbrace{\alpha_n^{(k)}p_n^{(t)}- c_n^{(t)}}_{\text{for }C}
    + 
    \underbrace{\beta_n^{(k)}y_n^{(t)} }_{\text{for }B}
    \le 0
    \\
    & 
    \underbrace{(\alpha')_n^{(k)}(p')_n^{(t)}
    - (c')_n^{(t)}}_{\text{for }C}
    + 
    \underbrace{(\beta')_n^{(k)}y_n^{(t)}}_{\text{for }B}
    \le 0 
\end{aligned}\tag{7, 8}
$$

Quick Start Binary Constraints

$$
\begin{aligned}
\text{For }B
\begin{cases}
    
    (x')^{(t)}_m + (z')^{(t)}_m
    \le 1
    \\[1.1em]
    (y')^{(t)}_m - (y')^{(t - 1)}_m - (x')^{(t)}_m + (z')_m^{(t)} 
    = 0
    \\[1.1em]
    
    \left(
        \sum_{\tau=t- T_{minu} + 1}^{t}
        (x')^{(\tau)}_m
    \right) - (y')^{(t)}_m \le 0 
    \\[1.1em]
    \left(
        \sum_{\tau = t - T_{mind} + 1}^{t}
        (z')^{(\tau)}_m
    \right)
    + (y')^{(t)}_m \le 1
\end{cases}
\end{aligned}\tag{9, 10, 11, 12}
$$

Capacity Constraints

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
    \underbrace{(sr)^{(t)}_n }_{\text{For} C}
    \underbrace{-R_n^{MT}y_n^{(t)}}_{\text{For }B} \le 0 
    \\&
    \underbrace{p_n^{(t)} + (sr)_n^{(t)} + (regu)_n^{(t)}}_{\text{for }C}
    \underbrace{- P^{\max}_ny_n^{(t)}}_{\text{for } B} \le 0
    \\&
    \underbrace{- p_n^{(t)} + (regd)^{(t)}_n}_{\text{for }C} 
    + \underbrace{y_n^{(t)}P_n^{min}}_{\text{for } B} \le 0
    \\
    & 
    \underbrace{(nsp)_n^{(t)}}_{\text{for }C} + \underbrace{\text{NSP}_n y_n^{(t)}}_{\text{for }B} \le \text{NSP}_n
    \\
    & 
    \underbrace{(regu)^{(t)}_n}_{\text{for }C}
    \underbrace{-\text{REGU}_ny_n^{(t)}}_{\text{for }B} 
    \le 0
    \\
    & 
    \underbrace{(regd)^{(t)}_n}_{\text{for }C} 
    \underbrace{-\text{REGD}_ny_n^{(t)}}_{\text{for }B} 
    \le 0
\end{aligned}\tag{13,..., 20}
$$

Quick start capacity constraints: Exactly the same as above, but wrt to the variable $p', x', y', z', sr', regd', regu', nsp'$. Whenever it's said "for C", it's still true, but "For B" is now "For G". 

Minimum Requirement Constraints


$$
\begin{aligned}
    \underbrace{- \sum_{n\in \mathcal G}^{}
    regu^{(t)}_n 
    -
    \sum_{m\in \mathcal G'}^{}
    (regu')^{(t)}_m }_{\text{for }C}
    \le 
    - \text{RREGU}^{(t)}
    \\
    \underbrace{- \sum_{n\in \mathcal G}^{}
    regd^{(t)}_n 
    -
    \sum_{m\in \mathcal G'}^{}
    (regd')^{(t)}_m }_{\text{for }C}
    \le 
    - \text{RREGD}^{(t)}
    \\
    \underbrace{- \sum_{n\in \mathcal G}^{}
    nsp^{(t)}_n 
    -
    \sum_{m\in \mathcal G'}^{}
    (nsp')^{(t)}_m }_{\text{for }C}
    \le 
    - \text{RNSP}^{(t)}
\end{aligned}\tag{30, \dots, 33}
$$

Battery Constraints: 


$$
\begin{aligned}
    \text{for C}
    \begin{cases}
        h_s^{(t + 1)} - h^{(t)}_s - \nu^+_s (g_s^+)^{(t)} + \nu^-_s(g_s^-)^{(t)} 
        \le 0 
        \\
        - h_s^{(t + 1)} + h^{(t)}_s + \nu^+_s (g_s^+)^{(t)} - \nu^-_s(g_s^-)^{(t)} 
        \le 0 
        \\
        h_s^{(t)} \le \overline{H}_s
        \\
        (g_s^+)^{(t)}  \le \overline{G}^+_s
        \\
        (g_s^-)^{(t)}  \le \overline{G}^-_s
    \end{cases}
\end{aligned}\tag{34, \dots, 38}
$$

Demand balance constraints

$$
\begin{aligned}
   \underbrace{ \sum_{n \in \mathcal G}^{}
    p_n^{(t)} + \sum_{m \in \mathcal G'}^{}(p')^{(t)}_m
    + \left(
        \sum_{s\in S}^{}(g^-_s)^{(t)} - (g^+_s)^{(t)}
    \right)}_{\text{for }C}
    -
    \underbrace{\sum_{b \in \mathcal B}^{}d_b^{(t)}}_{\text{For }F}
    &\le 0
    \\
    \underbrace{ -\sum_{n \in \mathcal G}^{}
    p_n^{(t)} - \sum_{m \in \mathcal G'}^{}(p')^{(t)}_m
    + \left(
        \sum_{s\in S}^{}(g^+_s)^{(t)} - (g^-_s)^{(t)} 
    \right)}_{\text{for }C}
    + 
    \underbrace{\sum_{b \in \mathcal B}^{}d_b^{(t)}}_{For F}
    &
    \le 0
    \\
    \underbrace{\sum_{b\in \mathcal B}^{}\sigma_b^{(l)}
    \left(
        \sum_{n \in G^b}^{}
            p_n^{(t)}
        + 
            \sum_{m \in (\mathcal G')^{b}}
                (p')_m^{(t)}
    \right) + 
    \sum_{s \in \mathcal S}^{}\mu_s^{(l)}
        \left(
            (g^-)_s^{(t)} - (g^+)_s^{(t)}
        \right)}_{\text{for }C}
    - 
    \underbrace{\sum_{b\in \mathcal B}^{}\sigma_b^{(l)}d_b^{(t)}}_{\text{For }F}
    &\le F^{(l)}
    \\
    \underbrace{- \sum_{b\in \mathcal B}^{}\sigma_b^{(l)}
    \left(
        \sum_{n \in G^b}^{}
            p_n^{(t)}
        + 
            \sum_{m \in (\mathcal G')^{b}}
                (p')_m^{(t)}
    \right)
    - \sum_{s \in \mathcal S}^{}\mu_s^{(l)}
        \left(
            (g^-)_s^{(t)} - (g^+)_s^{(t)}
        \right)}_{\text{for }C}
    + 
    \underbrace{\sum_{b\in \mathcal B}^{}\sigma_b^{(l)}d_b^{(t)}}_{\text{For }F}
    &\le F^{(l)} 
\end{aligned}
\tag{39, \dots, 44}
$$

The set $\mathcal G^b, (\mathcal G')^b$ are the primary and secondary generators on each of the busses. 


