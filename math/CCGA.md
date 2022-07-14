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


**FSP: Feasible Search Problem**
> * Given a demand $(d^{(k)})^*$ from FMP and a $\bar w$ from MP. 
> * It Controls: $u, v, q$. 
> * Objective: Minimizes the feasibility slacks variables, across all demands configurations. Searches for the lower bound for the feasibility slacks. 
> * Supports: 
>   * Produces the $(q^{(k)})^*, (u^{(k)})^*$ for the secondary decision variables. 
 

---
### **The Algorithm**

* Use MSP to identify a feasible system of variables: $q, u, w, \hat d$.
* Use FMP: 
  * Give: $M, \hat d$.  
  * Get: $\bar w, \bar \gamma$.
* Initiate $\epsilon$
* $k:=1$
* Initialize $q^{(0)}$ as random binary vector. 
* For Looping over some fixed amount of times: 
  * For $k$ incrementing one by one to a certain limit:  
    * use FMP
      * Give $q^{(k - 1)}, \bar w$ to FMP
      * Get $U$ the upper feasibility bound. 
      * Get $(\rho^+)^{(k)}, (\rho^-)^{(k)}, (d^{(k)})^*$
    * Use FSP: 
      * Give: $\bar w, (d^{(k)})^*$
      * Get lower bound $L$. 
    * If $L \ge 0$
      * Use MP
        * Introduce Cut using: $(u^{(k)}), (q^{(k)}), (\rho^+)^{(k)}, (\rho^-)^{(k)}$
        * Update $\bar w, \bar \gamma$
        * Breaks
    * elseif $U = 0$ OR $L = 0, U \le \epsilon$
      * Return (Terminates Everything)
    * Else: 
      * We don't know. 
  * Change the $\bar \gamma, \bar w$ in the FMP. 

**Remarks**

All the CCGA entities are kept when the algorithm is running over. They maintain the changes throughout the algorithm.  


