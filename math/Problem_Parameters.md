### **Intro**

There are 2 parts to it, Data and problem parameters. We need to convert the data to parameters appropriately. 

### **Global Constant**

$$
\begin{aligned}
    T_{end}, \mathcal B, \mathcal G, \mathcal G', \mathcal G^b, \mathcal G^{'b}, \mathcal L
\end{aligned}
$$

### **Generator Constants**

$$
\begin{aligned}
    NSP, R^{MT}
\end{aligned}
$$

---
### **Data Modeling**

General Format: 
- Data Entity
  - Attribute in code (Full name): Attribute in data sheet



- Generator: 
  - Index 
  - $P_{min}$: Pmin
  - $P_{max}$: Pmax
  - $T_{minu}$: Grow
  - $T_{mind}$: Hide
  - $regu$: 50% Pmax
  - $regd$: 50% Pmax
  - $sr$: 50% Pmax
  - $NSP$: 50% Pmax, CONST
  - $R^{MT}$: 50% Pmax, CONST
  - $\alpha^{(jb)}$: See Explainations
  - $\beta^{(jb)}$: See Explainations
  - $RREGU^t = 0$
  - $RREGD^t = 0$
  - $SD = 0$
  - $SU$: mu
  - $RU$: ????
  - $\overline{RU}$: Pmax
  - $\overline{RD}$: Pmax
  - $\overline{RD}$: ????
  - $H_n$: ????
  - $S_n$: ????


There are many instance of generators. $\alpha, \beta$ are related to the cost curve. "CONST" denotes that this it's the same for all generators. 

- Buss
  - $\hat{d} = 30$
  - $\bar{d}$



- Transmission Line
  - index
  - $F$: Limit
  - $SF$: Shear Factor
  - $X$

Transmission line goes between buses, they are like edges on the graph. 

---
### **Parameters: $\alpha_n^{jb}, \beta_n^{jb}$**

The set $jb$ are uniform discretizations of a domain of a quadratic function that models the cost for each of the bus, the domain is: $[P_{min}, P_{max}]$. $\alpha, \beta$ are slope and intercept for each of the interval on the domain. 



--- 
### **Parameters: $SF_l$**


