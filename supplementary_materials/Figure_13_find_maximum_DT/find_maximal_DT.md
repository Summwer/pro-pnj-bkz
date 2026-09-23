# Maximum of $\Delta_T(n,q)$

Suppose we have a LWE instance with parameters $n,q,m,\sigma$.

Assume $m = n $ and $\sigma =1$, then we have only two parameters $n$ and $q$. 

Since with $n$ (resp. $q$) growing, $\Delta_T(n,q) = T(S_0) - T(S_{op})$ first increase from 0, then decrease back to 0, there exists a pair $(n,q)$ maximizes $\Delta_T(n,q)$. 







### Design Experiment

Use cost model: MATZOV

$n = 100$~400, $q = 2^{10}$ ~ $2^{30}$.  $m = n$



补充测试：

$n = 500 \sim 600$, $q = 2^{20} \sim 2^{35}$. 



相同n下最优gap的q取值：

n = 450, $q = 2^{32}$, diff = 7.51



n = 500, $q=2^{35}$, diff = 8.41

n = 550, $q=2^{37}$, diff = 9.35

n = 600, $q=2^{41}$, diff = 10.34

diff≈0.01887n−0.98

q ≈ 2^{0.198 · n^{0.832}}

