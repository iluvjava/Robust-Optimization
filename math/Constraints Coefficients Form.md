### **Intro**
We seek to establish the standard matrix for the simple std form for the Secondary and the primary constraints groups. 





---- 
### **Notations**
Subscript Indexer: 

$$
[\mathcal G \times \mathcal H]:= \{(x, y) \in \mathcal G \times \mathcal H\}
$$


---
### **1st-Stage Constraints**
Relevant Matrix Coefficient Matrix: $A$. Define indexer and indexing var: $(t, n) \in [\mathcal G\times \mathcal T]$. 


### **First Stager Generator**

**Equality Constraints:**

$$
\begin{aligned}
    C3[(t, n)\in [\mathcal G\times \mathcal T], =0]:= 
    \begin{cases}
        y_n^{(t)}- y^{(t - 1)}_n
        \\
        -x^{(t)}_n
        \\
        + 
        z^{(t)}_n
    \end{cases}
\end{aligned}
$$

**Inequality Constraints:**

$$
\begin{aligned}
    C2[(t, n)\in [\mathcal G\times \mathcal T], \le 1]
    & :=
    \begin{cases}
        x^{(t)}_n
        \\
        z^{(t)}_n
    \end{cases}
    \\
    C4[(t, n) \in [\mathcal G\times \mathcal T],\le 0]
    &:= 
    \begin{cases}
        x^{(\tau)}_n & t - T_{n}^{(\text{minu})} + 1 \le \tau \le t
        \\
        - y^{(t)}_n & 
    \end{cases}
    \\
    C5[(t, n)\in[\mathcal G\times \mathcal T], \le 1]
    & := 
    \begin{cases}
        z^{(\tau)}_n & t - T_n^{(\text{mind})}+ 1\le \tau\le t
        \\
        -y^{(t)}_n & 
    \end{cases}
\end{aligned}
$$


---
### **2nd-Stage Constraints**

#### **Fuel Constraints**

**Variables**: 

$$
\begin{aligned}
    c_n^{t}
\end{aligned}
$$

Inequality Constraints
$$
\begin{aligned}
    C6[\empty, \le \phi] := 
    \begin{cases}
        \text{SU}_nx_n^{(t)} 
        &
        (t, n)\in [\mathcal G\times \mathcal T]
        \\
        \text{SD}_nz_n^{(t)} 
        &
        (t, n)\in [\mathcal G\times \mathcal T]
        \\
        c_n^{(t)}
    \end{cases}
\end{aligned}
$$


#### Quick Start Binary

#### Capacity Constraints

#### Quick Start Capacity Constraints

#### Demand Balance Constraints

#### Minimum Requirements Constraints

#### Intertia Constratins

---
### **Matrix A**

RHS: 
$$
\begin{aligned}
    w &= 
    \begin{bmatrix}
        x_n^{(t)} \\ y_n^{(t)} \\ z_n^{(t)}
    \end{bmatrix}
    \\
    & 
    \begin{cases}
        n \in \mathcal G
        \\
        t \in \mathcal T
    \end{cases}
    \\
    x,y,z
    &\in 
    \{0, 1\}^{|\mathcal G|\times |\mathcal T|}
\end{aligned}
$$

---
### **Matrix B**

RHS:

$$
\begin{aligned}
    w
\end{aligned}
$$


---
### **Matrix C**

RHS: 

$$
\begin{aligned}
    u 
    :=& 
    \begin{bmatrix}
        \begin{bmatrix}
             c_n^{(t)}
            \\
            (c')_m^{(t)}
            \\
            p_n^{(t)}
            \\
            (p')_m^{(t)}
            \\
            (sr)_n^{(t)}
            \\
            (sr')_m^{(t)}
            \\
            (regu)_n^{(t)}
            \\
            (regu')_m^{(t)}
            \\ 
            (nsp)^{(t)}_n
            \\
            (nsp')^{(t)}_m
            \\
            (regd)_n^{(t)}
            \\
            (regd')_m^{(t)}
        \end{bmatrix}
        \\
        \\
        \begin{bmatrix}
            (fft)_b^{(t)}
            \\
            (frr)_b^{(b)}
            \\ 
            (rfrr)_b^{(t)}
            \\
            (tffr)_b^{(t)}
            \\ 
            d^{(t)}_b
            \\
            (inx)^{(t)}_{sg}    
        \end{bmatrix}
    \end{bmatrix}
    \\
    & 
    \begin{cases}
        n \in \mathcal G
        \\
        m \in \mathcal G'
        \\
        t \in \mathcal T
        \\
        b\in \mathcal B 
        \\
        sg \in \mathcal N
    \end{cases}
\end{aligned}
$$

---
### **Matrix G**

RHS: 


$$
\begin{aligned}
    & 
    q := \begin{bmatrix}
        \begin{bmatrix}
            (x')_m^{(t)}
            \\
            (y')_m^{(t)}
            \\
            (z')_m^{(t)}
        \end{bmatrix}
        \\
        \begin{bmatrix}
            \delta_{sg}^{(t)}
        \end{bmatrix}
    \end{bmatrix}
    \\
    & 
    \begin{cases}
        m \in \mathcal G'
        \\
        t \in \mathcal T
        \\
        sg \in \mathcal N
    \end{cases}
\end{aligned}
$$