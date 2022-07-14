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
>   * Feasibility cut given: $u, q, d, \rho^+, \rho^-$ from FMP.
>   * Producing $\bar w, \bar \gamma$ as candicates primary decision variable and bounds around the average demands. 


**FSP: Feasible Search Problem**
> * Given a demand $(d^{(k)})^*$ that we wish to satisfy. 
> * It Controls the secondary variables, it tries to satisfies all the demands configurations. Secondary variables: $u, v, q$
> * Objective: Minimizes the feasibility slacks variables, across all demands configurations. Searches for the lower bound for the feasibility slacks. 
> * Supports: 
>   * Produces the $q^{(k)}, u^{(k)}$ for the secondary decision variables. 


**FMP: Upper Bound Searcher**
> * Given: $\hat d, w, q$
> * Controls over the secondary continuous and discrete variables: $u, v$ for all the demands profiles, and the demands itself for all profiles. 
> * It needs $q^{(k)}$ as an input for the the $k$ iteration of the CCGA. 
> * Objective: It tries to maximize the feasibility slack to locate a potential upper bound. 
> * Supports:
>   * Produces: $v^{(k)}, u^{(k)}, (\rho^+)^{(k)}, (\rho^-)^{(k)}$
>   * Accepts a fixed $q^{(k)}$, a specific secondary configurations. 



---
### **The Algorithm**
MSP: Identify candiate $\bar{w}, \bar{\gamma}$. 

**For** Any $q^{(k)}\in Q$ construct FMP:
  * FMP: Finds demand $(d^{(k + 1)})^*$ by maximizing slack variable $v$ for the secondary constraints. Demand and slack $v$ is a decision variable for the FMP. 
    * Update Upper bound: $UB$, upon first initialization, choose some random $Q$ regardless. 
  * FSP: Minimize the slacks variable for the secondary constraints and try to make it feasible given the demands $(d^{(k + 1)})^*$ found by the FMP. 
    * Update the Lower bound: $LB$
    * Obtain feasible secondary decision variables: $(u^{(k)})^*, (v^{(k)})^*, (q^{(k)})^*$. 
  * If $UB - LB \le \epsilon$, **terminates** CCGA. 
  * Create new decision variable: $u^{(k + 1)}, v^{(k + 1)}, \lambda^{(k + 1)}$ and add some constraints to the FMP. 
  * $k:= k + 1$
  * If, $UB - LB \approx 0$, then candidate $\bar{w}, \bar{\gamma}$ is feasible. else it has to be that $UB - LB > 0$, Perform feasibility cut to MP, and REPEAT the whole process, but keep the $q^{(k)}$ for the FMP, FSP in the CCGA. 



---
### **The Algorithm**

* Use MSP to identify a feasible system of variables: $q, u, w, \hat d$. 
* Use FMP: 
  * Give: $M, \hat d$.  
  * Get: $\bar w, \bar \gamma$. 
* Initiate $\epsilon$
* $k:=1$
* $q^{(k)} := q$
* Use FSP: 
  * Give: $\bar w, \hat d$
  * Get $(d^{(k)})^*$
  * Get lower bound $L$. 
* For Looping over some fixed amount of times: 
  * For $k$ incrementing one by one:  
    * Give $q^{(k)}$ to FMP
      * Get $U$ the upper feasibility bound. 
      * Get $(v^{(k)})^*, (u^{(k)})^*, (\rho^+)^{(k)}, (\rho^-)^{(k)}$
    * If $U - L \le \epsilon$
      * Breaks

**Remarks**

All the CCGA entities are kept when the algorithm is running over. They maintain the changes throughout the algorithm.  


