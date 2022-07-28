### **Master Problems Decision Variables and Constraints: Candidates Search**

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
### **FMP: Limited FCP; Feasibiilty Upper Bound, Adverserial Demands**

This is the objective of FMP: Fixing the secondary and primary decision variable, how possible is it to violate the constraints giving the freedom to change the demand within the bond specified by MP? In brief, the FMP algorithm is going to hold control for the following variables: 

* $v^{(k)}$, the dual slack variable for the secondary constraints paramaterized by CCGA iterations. 
* $d$, the demand vector. 
* $u$, constinous decision variable for the secondary problem. 
* $\lambda, \rho, \xi$: auxilary decision variables for colinear constraints from Complementary Slack of KKT.  



**Secondary Decision Variables**

Parameters: $u, q$, continuous and discrete variables for secondary problem. All discrete variable are binary. Because the variables inside are not "uniform" (different index set), the variable will be implemented by a Dictionary mapping from the name of the variable to a vector representing that variable. 

Continuous subvariables of $u$ are: 

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
\end{aligned}
$$

Discrete subvariables of $q$ are going to be:
$$
\begin{aligned}
    & \delta_{sg}^{(t)}, (x')_n^{(t)}, (y')_n^{(t)}, (z')_n^{(t)}
    \\
    & n\in \mathcal G', t\in \mathcal T, sg\in \mathcal N
\end{aligned}
$$

However, the CCGA creates another parameters for the $q$ indicating each step of the $q$ veraible in CCGA: $q^{(k)}$. And as a results, all inner variables are paramaterized too, giving us: 

$$
\begin{aligned}
    & \delta_{sg}^{(t, k)}, (x')_n^{(t, k)}, (y')_n^{(t, k)}, (z')_n^{(t, k)}
    \\
    & n\in \mathcal G', t\in \mathcal T, sg\in \mathcal N
\end{aligned}
$$

Where $k$ is related to the step of CCGA. 


**Dual Slack Variables and Adverserial Demands**

Variablse for the FCP for KKT conditions (Complementary Slack) are introduced. More specifically the variable $\xi$ and $\rho$. The FMP will optimize on $\eta$, an lower bound for slack, and then variables $\lambda, \xi, \rho, d$ are free variables. 



**Variables Modeling**

Both composite variables $u, q$ are going to be implemented via a dictionary mapping from the name of the variable as symbols for the key and maps to the instance for the variable in the code. 

**Constraints Matrix**
It's part of the representation of the constraints. 

* The Coefficient matrix for $u$ is denoted as $C$. 
* The Coefficient matrix for $q$ is denoted as $G$. 
* The vector $h$ is the second stage constraint RHS vector. 

The constraints matrix will be implemented implicitly as a set of paramatarized group of functions that perform the modifications for the FCP problem. Some tricks are applied to handle the uncertainty interval using the KKT conditions, giving us the format of the problem in the final form: 

$$
\begin{gathered}
    \psi(\bar{w}, \bar{\gamma}) = \eta\\
    \eta \leq \sum_{j=1}^{J} v_{j}^{r} r=1, \ldots, k \\
    C u^{r}-v^{r} \leq d+h-B \bar{w}-G q^{r} \quad r=1, \ldots, k \\
    \left(\lambda^{r}\right)^{T} C \leq 0, \quad r=1, \ldots, k \\
    \lambda^{r} \geq-1 \quad r=1, \ldots, k \\
    \lambda^{r}\left(d+h-B \bar{w}-G q^{r}\right)=\sum_{j=1}^{J} v_{j}^{r}, \quad r=1, \ldots, k \\
    \hat{d}-\bar{\gamma} \leq d \leq \hat{d}+\bar{\gamma}
\end{gathered}
$$

**Overall format**: 

$$
\begin{aligned}
    Bw + Cu + Gq \le d + h
\end{aligned}
$$

**Paramaterization of Constraints for FMP**


The constraints are paramaterized by candidates from the Master Problem: $\bar{w}, \bar{\gamma}, k$, and constants $q\in Q$. $k$ is the current iterations of the CCGA algorithm, and it can be accumlated for the FMP algorithm. 


---
### **FSP: Feasibility Lower Bound, Recovering Feasible Parameters for Adverserial Demands**

Given demands output from FSP, $\bar{w}, \bar{\gamma}$ from MP, FSP determines the feasible solutions by changing secondary decision variable $u, q$, and optimizing for the slack variable $v$. It's asking the question that, given the demands is it possible to find the most feasible solution it by choosing the approriate secondary variables? 







---
### **FMP Interfacing with Others**

**Integration of CCGA**

> To integrate the decision variable with CCGA, we just put all the decision variables during the iterations of the CCGA into a list, list index is the same as the superscript for the variables: $v^{(k)}$


