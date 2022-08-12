### **Introduction**

The strong duality constraints presented us with the following equality constraints for each $k$ iterations of the CCGA inner loop: 

$$
\begin{aligned}
   (\lambda^{(k)})^T(Hd +h - Gq^{(k)*} - B\bar{w}) = \sum_{j = 1}^{J}v_j^{(k)}
\end{aligned}
$$

There are 2 ways we are trying to taming the constraints: 
* Formulating it using the corner point assumptions 
* Formulating it using convex relaxation for non-convex constraints: McCormick Envelope. 


---
### **Bi-Linear Reformulations via Corner Point Assumptions**


We assume that the demands config that achieves the optimal for the primals are the the corner points of the hyper cube: $\bigcup_{b\in \mathcal{B}}[\hat d_b - \bar\gamma, \hat d_b + \bar\gamma]$, to achive that we consider demands decision variable in the form of: 

$$
\begin{aligned}
    & d := \bar\gamma \mathbf \rho^+ + \hat d - \bar{\gamma}\rho^-
    \\
    & \rho^+, \rho^- \in \{0, 1\}^{|\mathcal B|}
\end{aligned}
$$
The decision variable $\rho^+, \rho^-\in \{1, 0\}^{|\mathcal B|}$. More importantly, they are not indexed by: $k$. We reconsider the cross term $(\lambda^{(k)})^THd$ which should give us: 

$$
\begin{aligned}
    &= (\lambda^{(k)})^THd
    \\
    &= (\lambda^{(k)})^TH\hat{d} + \bar\gamma(\lambda^{(k)})^TH\rho^{+} - \bar\gamma(\lambda^{(k)})^TH\rho^{+}
\end{aligned}
$$

Given all choices of $\rho^+, \rho^-$, the equality should remains true by assumptions. Then we consider another matrix of decision variable $\Xi^+_k, \Xi^-_k$ such that $\sum_{(i, j)\in [J]\times [B]}H\circ \Xi^+_k = (\lambda^{(k)})^TH\rho^+$. Then we have: 

$$
(\lambda^{(k)})^THd = 
(\lambda^{(k)})^T H\hat d + \bar\gamma \sum_{(i, j)\in [J]\times [B]}^{}(H\circ (\Xi^+_k - \Xi^-_k))_{i, j}
$$

We use $(\Xi^+_k)_{j, b}$ to denote a particular decision variables indexed by $(j, b)$ inside of the matrix. The constraints that enfoce such equality for all corner points demands is formulated as: 

$$
\begin{aligned}
    \forall (j, b) :H_{j, b}\neq 0, \forall k: &
    \begin{cases}
        \lambda_j^{(k)}  - (1 - \rho^+_{b}) 
        \le (\Xi_k^+)_{j, b} 
        \le \lambda^{(k)}_j + (1 - \rho^+_{b})
        \\
        -\rho^+_b \le (\Xi_k^{+})_{j, b} \le \rho^+_b
        \\
        \lambda_j^{(k)}  - (1 - \rho^-_{b}) 
        \le (\Xi_k^-)_{j, b} 
        \le \lambda^{(k)}_j + (1 - \rho^-_{b})
        \\
        -\rho^-_b \le (\Xi_k^{-})_{j, b} \le \rho^-_b
    \end{cases}
    \\
    & \rho^+_b + \rho^-_b = 1 \quad \forall b
    \\
    & \rho^+, \rho^- \in \{0, 1\}^{|\mathcal B|}
\end{aligned}
$$


---
### **Addressing the Uneven Demands Issues**



---
### **McCormick Envelope: Convex Relaxations to Non-Convex Constraint**

The over leaf project link is [here](https://www.overleaf.com/project/62ed834707a735792c56b4ab). 