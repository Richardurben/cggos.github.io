---
layout: post
title:  "矩阵分解与线性方程组"
date:   2019-01-28
categories: Maths
tags: [Matrix, Least-squares Minimization]
---

[TOC]

# Matrix

## 正规矩阵

$$
A^H A = A A^H
$$

* 酉矩阵 $A^H A = A A^H = I$
  * 正交矩阵 $A^T A = AA^T = I$
    * 特殊正交矩阵：$\det(A) = +1$，是一个旋转矩阵
      * Givens矩阵：$G(i,j,\theta)$
    * 瑕旋转矩阵：$\det(A) = -1$，瑕旋转是旋转加上镜射
      * Householder矩阵：$H = I - 2 uu^T$

* 厄米特矩阵 $A^H = A$
  * 实对称矩阵 $A^T = A$

* 反厄米特矩阵 $A^H = -A$
  * 实反对称矩阵 $A^T = -A$


## 正定与半正定矩阵

在线性代数里，正定矩阵(positive definite matrix)有时简称为 **正定阵**。

对任意的非零向量 $\boldsymbol{x}$ 恒有 **二次型**  

$$
f = \boldsymbol{x}^T \boldsymbol{A} \boldsymbol{x} > 0
$$

则称 $f$ 为 **正定二次型**，对应的矩阵为 **正定矩阵**；若 $f \ge 0$，则 对应的矩阵为 **半正定矩阵**。


### 直观理解

令 $Y=MX$，则 $X^T Y > 0$，所以  

$$
cos(\theta) = \frac{X^T Y}{\|X\|\|Y\|} > 0
$$

因此，从另一个角度，正定矩阵代表一个向量经过它的变化后的向量与其本身的夹角 **小于90度**

### 判别对称矩阵A的正定性

* 求出A的所有特征值
  * 若A的特征值均为正数，则A是正定的；若A的特征值均为负数，则A为负定的。
* 计算A的各阶顺序主子式
  * 若A的各阶顺序主子式均大于零，则A是正定的；若A的各阶顺序主子式中，奇数阶主子式为负，偶数阶为正，则A为负定的。


## 其他

* **对角阵**：任意正规矩阵 都可以经过 正交变换 变成 对角矩阵
* **上（下）三角阵**
* **可逆矩阵（非奇异矩阵）**


## 矩阵变换

### 正交变换

### 吉文斯旋转（Givens rotation）

在数值线性代数中，吉文斯旋转（Givens rotation）是在两个坐标轴所展开的平面中的旋转。

吉文斯旋转 矩阵：
$$
G(i, j, \theta)=\left[\begin{array}{ccccccc}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & -s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{array}\right]
$$

当一个吉文斯旋转矩阵 G(i,j,θ)从左侧乘另一个矩阵 A 的时候，GA 只作用于 A 的第 i 和 j 行。

$$
\left[\begin{array}{cc}
c & -s \\
s & c
\end{array}\right]\left[\begin{array}{l}
a \\
b
\end{array}\right]=\left[\begin{array}{l}
r \\
0
\end{array}\right]
$$

明确计算出 θ 是没有必要的。我们转而直接获取 c, s 和 r。一个明显的解是

$$
\begin{aligned}
&r \leftarrow \sqrt{a^{2}+b^{2}}\\
&c \leftarrow a / r\\
&s \leftarrow-b / r
\end{aligned}
$$

Givens 旋转在数值线性代数中主要的用途是在向量或矩阵中介入零。例如，这种效果可用于计算矩阵的 QR分解。超过Householder变换的一个好处是它们可以轻易的并行化，另一个好处是对于非常稀疏的矩阵计算量更小。

### Householder 变换

豪斯霍尔德变换（Householder transformation）或译“豪斯霍德转换”，又称初等反射（Elementary reflection），这一变换将一个向量变换为由一个超平面反射的镜像，是一种线性变换。其变换矩阵被称作豪斯霍尔德矩阵，在一般内积空间中的类比被称作豪斯霍尔德算子。超平面的法向量被称作豪斯霍尔德向量。

<div align=center>
  <img src="../images/maths/householder_reflection.jpg">
</div>

豪斯霍尔德变换可以将向量的某些元素置零，同时保持该向量的范数不变。例如，将非零列向量 $\mathbf{x}=\left[x_{1}, \dots, x_{n}\right]^{T}$ 变换为单位基向量 $\mathbf{e}=[1,0, \ldots, 0]^{T}$ 的豪斯霍尔德矩阵为

$$
\mathbf{H}=\mathbf{I}-\frac{2}{\langle\mathbf{v}, \mathbf{v}\rangle} \mathbf{v} \mathbf{v}^{H}
$$

其中豪斯霍尔德向量 $\mathbf{v}$ 满足：

$$
\mathbf{v}=\mathbf{x}+\operatorname{sgn}\left(\mathbf{x}_{1}\right)\|\mathbf{x}\|_{2} \mathbf{e}_{1}
$$

# Matrix Decomposition

## EVD (Eigen Decomposition)

$$
A = VDV^{-1}
$$

* $A$ 是 **方阵**；$D$ 是 **对角阵**，其 **特征值从大到小排列**；$V$ 的列向量为 **特征向量**
* 若 $A$ 为 **对称阵**，则 $V$ 为 **正交矩阵**，其列向量为 $A$ 的 **单位正交特征向量**

## SVD (Singular Value Decomposition)

<div align=center>
  <img src="../images/maths/svd.jpg">
</div>

$$
A = UDV^T
$$

* $A$ 为 $m \times n$ 的矩阵；$D$ 是 **非负对角阵**，其 **奇异值从大到小排列**；$U$、$V$ 均为 **正交矩阵**

SVD分解十分强大且适用，因为任意一个矩阵都可以实现SVD分解，而特征值分解只能应用于方阵。

### 奇异值与特征值

$$
AV = UD \Rightarrow Av_i = \sigma_{i} u_i \Rightarrow
\sigma_{i} = \frac{Av_i}{u_i} \\[4ex]
A^T A = (V D^T U^T) (U D V^T) = V D^2 V^T \\[2ex]
A A^T = (U D V^T) (V D^T U^T) = U D^2 U^T
$$

* $A^T A$ 的 **特征向量** 组成的是SVD中的 $V$ 矩阵
* $A A^T$ 的 **特征向量** 组成的是SVD中的 $U$ 矩阵
* $A^T A$ 或 $A A^T$ 的 **特征值** $\lambda_i$ 与 $A$ 的 **奇异值** $\sigma_i$ 满足 $\sigma_i = \sqrt{\lambda_i}$

### PCA

**奇异值减少得特别快**，在很多情况下，**前10%甚至1%的奇异值的和就占了全部的奇异值之和的99%以上**，可以用 **最大的 $k$ 个的奇异值和对应的左右奇异向量** 来近似描述矩阵  

$$
A_{m \times n} = U_{m \times m} D_{m \times n} V_{n \times n}^T
\approx U_{m \times k} D_{k \times k} V_{k \times n}^T
$$


## LU Decomposition

$$
A = LU
$$

* $A$ 是 **方阵**；$L$ 是 **下三角矩阵**；$U$ 是 **上三角矩阵**

### PLU 分解

$$
A = PLU
$$

事实上，PLU 分解有很高的数值稳定性，因此实用上是很好用的工具。

### LDU 分解

$$
A = LDU
$$

## Cholesky Decomposition

* [Cholesky decomposition](https://rosettacode.org/wiki/Cholesky_decomposition)

$$
A = LDL^T
$$

* $A$ 是 **方阵**，**正定矩阵**；$L$ 是 **下三角矩阵**

classic:

$$
A = LL^T \\[2ex]
A^{-1} = (L^T)^{-1} L^{-1} = (L^{-1})^T L^{-1}
$$

### LDLT

* [LDLT for Checking Positive Semidefinite Matrix's Singularity](http://simbaforrest.github.io/blog/2016/03/25/LDLT-for-checking-positive-semidefinite-matrix-singularity.html)

## QR Decomposition

* [QR decomposition](https://rosettacode.org/wiki/QR_decomposition)

$$
A = QR
$$

* $A$ 为 $m \times n$ 的矩阵；$Q$ 为 **酉矩阵**；$R$ 是 **上三角矩阵**


# 线性方程组

<div align=center>
  <img src="../images/maths/matrix_function_solve.jpg">
</div>

## 非齐次线性方程组

$$
A_{m \times n} x = b_{m \times 1}
$$

在非齐次方程组中，A到底有没有解析解，可以由增广矩阵来判断：

* r(A) > r(A | b) 不可能，因为增广矩阵的秩大于等于系数矩阵的秩
* r(A) < r(A | b) 方程组无解；
* r(A) = r(A | b) = n，方程组有唯一解；
* r(A) = r(A | b) < n，方程组无穷解；

### 非齐次线性方程组的最小二乘问题

$$
x^* = \arg \min_x{\|Ax - b\|}_2^2
$$

m个方程求解n个未知数，有三种情况：

* m=n
  - 且A为非奇异，则有唯一解 $x=A^{-1}b$
* m>n，**超定问题（overdetermined）**
  - $x=A^{+}b$
* m<n，**欠定问题（underdetermined）**


## 齐次线性方程组

$$
A_{m \times n} x = 0
$$

**齐次线性方程 解空间维数=n-r(A)**

* r(A) = n
  - A 是方阵，该方程组有唯一的零解
  - A 不是方阵(m>n)，解空间只含有零向量
* r(A) < n
  - 该齐次线性方程组有非零解，而且不唯一（自由度为n-r(A))


### 齐次线性方程组的最小二乘问题

$$
\min{\|Ax\|}_2^2 \quad s.t. \quad \|x\| = 1
$$

* 最小二乘解为 **矩阵 $A^TA$ 最小特征值所对应的特征向量**
* $\text{EVD}(A^{T}A)=[V,D]$，找最小特征值对应的V中的特征向量
* $\text{SVD}(A)=[U,S,V]$，找S中最小奇异值对应的V的右奇异向量


# with Eigen

* [Eigen 3.2稀疏矩阵入门](https://my.oschina.net/cvnote/blog/166980)

* [Catalogue of dense decompositions](https://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html)

* [Benchmark of dense decompositions](https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html)

* [Solving linear least squares systems](https://eigen.tuxfamily.org/dox/group__LeastSquares.html)
  - SVD decomposition (the most accurate but the slowest)
  - QR decomposition
  - normal equations (the fastest but least accurate)

* [Solving Sparse Linear Systems](http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html)
