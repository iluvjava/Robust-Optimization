### **CCGA Entities**


**MP: Main Problem**
> * Given: N/A
> * Controls over: $d, u, q, w, \gamma$. 
> * Objective: Varies
> * Purpose: Check if possible solutions exists and search for reasonable $\hat d$. 
> Supports:
>   * Produces: $d, u, q, w$ that is feasible, with some random objectives. 


**MSP: Master Problem:**
> * Given upper bound $M$ for $\gamma$. $\hat d$ that is at least feasible for the system. 
> * Purpose: set up the demands intervals and the primary binary decision variables. 
> * Controls $w, \gamma$
> * Maximizes $\gamma$
> * Supports: 
>   * Feasibility cut given: $u, q, \rho^+, \rho^-$ from FMP.
>   * Producing $\bar w, \bar \gamma$ as candicates primary decision variable and bounds around the average demands. 


**FMP: Upper Bound Searcher**
> * Given: $\hat d, \bar w, \bar \gamma$, and a $q$ that is arbitrary at the first time, later on it needs to be something be from FSP, previus iteration. 
> * Controls: $u, v$ for all the demands profiles, and the demands itself for all profiles. 
> * It needs $q^{(k)}$ as an input for the the $k$ iteration of the CCGA. 
> * Objective: It tries to maximize the feasibility slack to locate a potential upper bound. 
> * Supports:
>   * Produces: $u^{(k)}, q^{(k)}, (\rho^+)^{(k)}, (\rho^-)^{(k)}$ for MP cut ONLY. 
>   * Produces: $(d^{(k)})^*$
>   * Accepts a fixed $q^{(k)}$, a specific secondary configurations. 
>   * Automatically make constraints for newly introduced $q$ at each iterations. 


**FSP: Feasible Search Problem**
> * Given a demand $(d^{(k)})^*$ from FMP and a $\bar w$ from MP. 
> * It Controls: $u, v, q$. 
> * Objective: Minimizes the feasibility slacks variables, across all demands configurations. Searches for the lower bound for the feasibility slacks. 
> * Supports: 
>   * Produces the $(q^{(k)})^*, (u^{(k)})^*$ for the secondary decision variables. 
 

---
### **The Algorithm**

* Use MP to identify a feasible system of variables: $q, u, w, \hat d$.


* Initiate $\epsilon$ 
* For Looping over some fixed amount of times:
  * Use MSP: 
    * Give: $M, \hat d$.  
    * Get: $\bar w, \bar \gamma$.
  * Initialize $(q^{(0)})^*$ as random binary vector. 
  * Create FMP with $(q^{(0)})^*, \bar w, \bar \gamma$. 
  * Create FSP with a $\bar w$
  * For $k = 1$, incrementing one by one to a certain limit:  
      * use FMP
        * Give $q^{(0)}$ if $k = 1$ else $(q^{(k - 1)})^*$
        * Get $U$ the upper feasibility bound. 
        * Get $(\rho^+)^{(k)}_*, (\rho^-)^{(k)}_*, (d^{(k)})^*$
    * Use FSP: 
      * Give: $(d^{(k)})^*$
      * Get lower bound $L$. 
      * Get $(u^{(k)})^*, (q^{(k)})^*$. 
    * If $L > 0$
      * Use MP
        * Introduce Cut using: $(u^{(k)})^*, (q^{(k)})^*, (\rho^+)^{(k)}_*, (\rho^-)^{(k)}_*$
        * Breaks out the inner forloop. 
    * elseif $U \le \epsilon$
      * Save the results somewhere. 
      * Return (Terminates Everything)
    * Else: 
      * Nothing to be done here. 

**Remarks**

All the CCGA entities are kept when the algorithm is running over. They maintain the changes throughout the algorithm.  


---
### **The Problem of Performing Cut on the Master Problem**

Originally cut perform for the master problem involves $\rho^\pm, u, q, \hat{d}$. The solution from FSP is the simply the following:  

$$
\begin{aligned}
	Cu^* + H(d^{(k + 1)})^* - v^*
\end{aligned}
$$


---
### **The Feasibility Problem**

To test whether there exists a feasible solution for the polytope $Ax \le b$, we consider the objective of the following optimization problem: 

$$
\begin{aligned}
	\min_{v\ge 0}\{ v:Ax - v\mathbf{1} \le b\}
\end{aligned}
$$

Where, we are only optimizing $v$ which si a scalar, and if it's possible to set $v$ into $0$. Here, $v$ doesn't need to be a vector. But when I implemented the CCGA, the main problem, for FMP, FSP are all treating $v$ as a vector. It's still possible to minimize wrt to a vector $v$ that has the same dimension as the constraints


---
### **Another Unatural Things That I Found Out**

Suppose that $\gamma = 0, (d^{(k + 1)})^* = (u^{(k + 1)})^* = \mathbf 0$, then the cut we introduced in the original paper(with $-H$ here being the $H$ in the original paper): 

$$
	C(u^{k+1})^*
		+ 
	G(q^{k+1})^*
		+ 
	H\big(\hat{d}+(\rho_*^+)^{k+1}\gamma\bm{1} - (\rho_*^-)^{k+1}\gamma\bm{1}\big)
		+ 
	Bw
		\leq
	h
$$

Should simplify into: 

$$
Bw  + H\hat d\le h
$$

And intuitively for this particular case, I would expect that there is a solution using just the Primary Generator such that it can satisifes the average demands. However, when I experimented with such a cut, the system seems to be infeasible, which is very weird. 


