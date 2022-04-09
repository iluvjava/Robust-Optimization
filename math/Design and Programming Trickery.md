### **Programming Tricks and Design Philosophies**

**The master constraints**: 

> Here we make use of the master constraints which simultaneous appears in multiple context of the algorithm, but paramaterized different. The constraint is: 
> $$
> \begin{aligned}
>   Bw + Cu + Gq + v &\le d + h  \quad \leftarrow \text{MP}
>   \\
>   B\bar{w} + Cu^{(r)} + Gq^{(r)} + v^{(r)} &\le d + h  \quad \leftarrow \text{FCP}
> \end{aligned}
> $$

* The MP controls: $w$; Others are fixed. 
* The FCP controls: $u, v, q$. 

