### **Intro**

We need to explicitly form the constraint matrix when programming it up. 

Given the variable name, such as $c, p, rfrr$, with superscript and subscript, we need to translate these information to the column position for its coefficient. In this we, we can construct a sparse matrices for the problem constraints. The ordering must be kept **consistent**. 


---
### **Primary Decision Variables**


$$
w = \begin{bmatrix}
    \\
    x \in \mathbb R^{|G|\times |\mathcal T|}
    \\
    y \in \mathbb R^{|G|\times |\mathcal T|}
    \\
    z \in \mathbb R^{|G|\times |\mathcal T|}
    \\[1.1em]
\end{bmatrix}, \gamma \in \mathbb R_+
$$

The variable $x, y, z$ has underscript $x_n^{(t)},y_n^{(t)},z_n^{(t)}$

$$
\begin{aligned}
    & (t, n, s) \mapsto n |\mathcal G| + t + s|\mathcal G||\mathcal T|
    \\
    & x_{(t, n)} = w_{(t, n)}
    \\
    & y_{(t, n)} = w_{(t, n) + |\mathcal G||\mathcal T|}
    \\
    & y_{(t, n)} = w_{(t, n) + |\mathcal G||\mathcal T|}
\end{aligned}
$$

Relevant Constraints Matrices: $A, B$, they models constraints: 

* $A$: 2 -> 5



---
### **Secondary Variables: Discrete**

The secondary discrete variables are all packed into the same vector $q^{(k)}$. But before that, we make a list of then so we know what they are. 

* Constraint matrix: $G$

* Compositve variable: $q$

$$
w' = \begin{bmatrix}
    \\
    x' \in \mathbb R^{|G'|\times |\mathcal T|}
    \\
    y' \in \mathbb R^{|G'|\times |\mathcal T|}
    \\
    z' \in \mathbb R^{|G'|\times |\mathcal T|}
    \\[1.1em]
\end{bmatrix}
$$

$$
\begin{aligned}
    & \delta_{sg}^{(t)}, (x')_n^{(t)}, (y')_n^{(t)}, (z')_n^{(t)}
    \\
    & n\in \mathcal G', t\in \mathcal T, sg\in \mathcal N
\end{aligned}
$$



