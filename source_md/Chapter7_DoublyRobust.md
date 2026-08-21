# ECON5280 Chapter 7 Doubly Robust Learning

<font size="5">Junlong Feng</font>

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/junlong-feng/econ5280/main?filepath=Lecture7_DR.ipynb)

## Outline

* Motivation: Can we estimate ATE under unconfoundedness? If yes, what's the best way?
* OLS no longer works: ATE can be consistently estimated by OLS under complete randomization but that's no longer true under unconfoundedness in general.
* Two ways to identify ATE: Expectation of CATE and inverse propensity score weighting.
* Doubly Robust Learning: AIPW, the best way to estimate ATE under unconfoundedness.
* Implementation: Everything can be done in R by a simple command.

## 1. ATE under Unconfoundedness and the Failure of OLS

In Chapter 5, we said if we have complete randomization, i.e.,
$$
D_{i}\perp (Y_{i}(1),Y_{i}(0),W_{i}),
$$
we can estimate ATE of $D\in \{0,1\}$ on $Y$ by linearly regressing $Y_{i}$ on $(1,D_{i})$ or $(1,D_{i},W_{i})$, regardless of the true model.

In Chapter 6, we learned a relaxed assumption which may entertain a wider range of data. This weaker assumption is unconfoundedness:
$$
D_{i}\perp (Y_{i}(1),Y_{i}(0))|W_{i}.
$$
In this chapter, we consider the following questions:

- Is ATE identified under unconfoundedness?
- If yes, can we still estimate it using OLS?
- If yes to the first question but no to the second, then how should we estimate ATE in this case?

The answer to the first question is obvious: Yes because we learned in Chapter 6 that under unconfoundedness, CATE is identified, and since ATE is the expectation of CATE, we can identify ATE as:
$$
ATE\overset{unconfoundedness}{=}\mathbb{E}[\mathbb{E}(Y_{i}|D_{i}=1,W_{i})-\mathbb{E}(Y_{i}|D_{i}=0,W_{i})].
$$
But how about the second question? Is regressing $Y$ onto $(1,D)$ or $(1,D,W)$ still a good idea? 

```R
library(lmtest)
library(sandwich)
set.seed(5280)
## DGP
n=5000
w=rnorm(n,0,1)
v=rnorm(n,0,1)
d=((w+v)>0) ## So the treatment satisfies unconfoundedness but not complete randomization
e=rnorm(n,0,1)
y=1+d+w+d*w^2+sin(w)+e  ## So, ATE=2

## OLS
model1=lm(y~d)
model2=lm(y~d+w)
coeftest(model1,vcov.=vcovHC(model1,type="HC0"))
coeftest(model2,vcov.=vcovHC(model2,type="HC0"))
```

From the simulation example, we can see that the estimates are not close to the true ATE. The reason is as follows. Consider the regression 
$$
Y_{i}=\beta_{0}+\beta_{1}D_{i}+\beta_{2}W_{i}+u_{i}.
$$
By the true model, $u_{i}=D_{i}W_{i}^{2}+sin(W_{i})+e_{i}$. By construction $D_{i}$ and $W_{i}$ are correlated, so $D_{i}$ and $u_{i}$ are correlated, leading to the failure of OLS/MM. However, we can verify that unconfoundedness actually holds in this case. 

To summarize: Under unconfoundedness, ATE is identified, so it is indeed possible to unbiasedly and consistently estimated it. Now the problem is that OLS/MM is no longer the proper estimator. We need to seek a different estimation approach. If we can find one, then we achieve something huge because ATE will not only be estimable under complete randomization: It can also be well-estimated under conditional randomization or even in observational data!

## 2. Two Identification Equations for ATE

It turns out that under unconfoudedness, ATE is not only identified by equation (3). We can identify it from a different angle.

Since $D$ is binary, if $W$ and $D$ are correlated (this is the only interesting case because otherwise, we have complete randomization and thus don't need to worry about estimation), then $W$ only affects $D$ via the following conditional probability:
$$
p(W_{i})\equiv\Pr(D_{i}=1|W_{i}).
$$
This function $p(\cdot)$ is called the **propensity score** of $D$. It tells us how $W_{i}$ affects the take-up probability of the treatment. In the special case of complete randomization, $p(W_{i})$ is a constant, not changing at all no matter what value $W_{i}$ takes.

Since the propensity score is the only possible chanel for $W$ to affect $D$, if $D$ is independent of $(Y(1),Y(0))$ conditional on $W$, then one should expect to see that $D$ is independent of $(Y(1),Y(0))$ conditional on $p(W)$ as well. Indeed, this is a theorem and can be formally proved (Rosenbaum and Rubin, 1983).

**Theorem**. If $D_{i}\perp (Y_{i}(1),Y_{i}(0))|W_{i}$ and $D_{i}\in\{0,1\}$, then $D_{i}\perp (Y_{i}(1),Y_{i}(0))|p(W_{i})$.

*Proof*. By the definition of independence, we are done if we can prove $\Pr(D_{i}=1|Y_{i}(0),Y_{i}(1),p(W_{i}))=\Pr(D_{i}=1|p(W_{i}))$ . Now,
$$
\begin{align*}
&\Pr(D_{i}=1|Y_{i}(0),Y_{i}(1),p(W_{i}))\\
=&\mathbb{E}(D_{i}|Y_{i}(0),Y_{i}(1),p(W_{i}))\\
\overset{tower\ property}{=}&\mathbb{E}[\mathbb{E}(D_{i}|Y_{i}(0),Y_{i}(1),W_{i})|Y_{i}(0),Y_{i}(1),p(W_{i})]\\
\overset{unconfoundedness}{=}&\mathbb{E}[\mathbb{E}(D_{i}|W_{i})|Y_{i}(0),Y_{i}(1),p(W_{i})]\\
=&\mathbb{E}(p(W_{i})|Y_{i}(0),Y_{i}(1),p(W_{i}))\\
=&p(W_{i}),
\end{align*}
$$
where the penultimate equality is by $\mathbb{E}(D_{i}|W_{i})=\Pr(D_{i}=1|W_{i})\equiv p(W_{i})$.

Meanwhile, 
$$
\Pr(D_{i}=1|p(W_{i}))=\mathbb{E}(D_{i}|p(W_{i}))=\mathbb{E}[\mathbb{E}(D_{i}|W_{i})|p(W_{i})]=p(W_{i}),
$$
where the second equality is by the tower property. Done.

The beauty of this theorem is that, instead of controlling the whole vector $W_{i}$, which is multi-dimensional in general, one can simply condition on the propensity score, which is a scalar.

Now we're ready to derive an alternative identification equation for ATE. The derivation of this new equation does not use the theorem above, but the form of this identification equation is a perfect example that knowing the propensity score is as good as knowing the whole $W_{i}$ vector.

**Theorem**. Under unconfoundedness, if $p(W_{i})\in (0,1)$ almost surely (which is called **overlap**), then
$$
ATE=\mathbb{E}\left(\frac{D_{i}Y_{i}}{p(W_{i})}-\frac{(1-D_{i})Y_{i}}{1-p(W_{i})}\right).
$$
*Proof*. Note that $D_{i}Y_{i}=Y_{i}(1)$ if $D_{i}=1$ and $0$ if $D_{i}=0$. Therefore, $D_{i}Y_{i}=D_{i}Y_{i}(1)$. Similarly, $(1-D_{i})Y_{i}=(1-D_{i})Y_{i}(0)$ . Substitute them into the right hand side:
$$
\begin{align*}
&\mathbb{E}\left(\frac{D_{i}Y_{i}}{p(W_{i})}-\frac{(1-D_{i})Y_{i}}{1-p(W_{i})}\right)\\
=&\mathbb{E}\left(\frac{D_{i}Y_{i}(1)}{p(W_{i})}\right)-\mathbb{E}\left(\frac{(1-D_{i})Y_{i}(0)}{1-p(W_{i})}\right)\\
\overset{LIE}{=}&\mathbb{E}\left(\frac{1}{p(W_{i})}\mathbb{E}\left(D_{i}Y_{i}(1)|W_{i}\right)\right)-\mathbb{E}\left(\frac{1}{1-p(W_{i})}\mathbb{E}\left((1-D_{i})Y_{i}(0)|W_{i}\right)\right)\\
\overset{unconfoundedness}{=}&\mathbb{E}\left(\frac{1}{p(W_{i})}\mathbb{E}\left(D_{i}|W_{i}\right)\mathbb{E}\left(Y_{i}(1)|W_{i}\right)\right)-\mathbb{E}\left(\frac{1}{1-p(W_{i})}\mathbb{E}\left((1-D_{i})|W_{i}\right)\mathbb{E}\left(Y_{i}(0)|W_{i}\right)\right)\\
=&\mathbb{E}\left(\mathbb{E}\left(Y_{i}(1)|W_{i}\right)\right)-\mathbb{E}\left(\mathbb{E}\left(Y_{i}(0)|W_{i}\right)\right)\\
\overset{by\ eq.(3)}{=}&ATE.
\end{align*}
$$

So now, we can write ATE either as equation (3) or equation (7). This motivates us to estimate ATE by jointly using both.

## 3 Augmented Inverse Propensity Score Weighting (AIPW) and Double Robustness 

### 3.1 An Oracle Estimator to Illustrate the Idea

For simplicity, let $f_{d}(W_{i})\equiv\mathbb{E}(Y_{i}|D_{i}=d,W_{i})$ . Then equation (3) says
$$
\tau\equiv ATE=\mathbb{E}(f_{1}(W_{i})-f_{0}(W_{i})).
$$
So it is natural to estimate ATE by
$$
\hat{\tau}^{*}_{ECATE}\equiv\frac{1}{n}\sum_{i=1}^{n}\left(f_{1}(W_{i})-f_{0}(W_{i})\right).
$$
Of course, in practice, we would never observe $f_{1}$ and $f_{0}$. So $\tau^{*}_{ECATE}$ is called an "oracle" estimator, an infeasible estimator where some true functions/parameters are involved. Oracle estimators are often helpful benchmarks to fix ideas and to focus on the main insight.

Meanwhile, by equation (7), another oracle estimator of ATE could be
$$
\hat{\tau}^{*}_{IPW}\equiv \frac{1}{n}\sum_{i=1}^{n}\frac{Y_{i}D_{i}}{p(W_{i})}-\frac{1}{n}\sum_{i=1}^{n}\frac{Y_{i}(1-D_{i})}{1-p(W_{i})}.
$$
$\hat{\tau}^{*}_{IPW}$ is also an oracle estimator because in most data, we do not known the true propensity score function $p$.

Now, let's think about how to combine (9) and (10). We propose the following *infeasible augmented inverse propensity score weighted* estimator:
$$
\begin{align*}
\hat{\tau}^{*}_{AIPW}\equiv& \frac{1}{n}\sum_{i=1}^{n}\left(f_{1}(W_{i})-f_{0}(W_{i})\right)\\
&+\frac{1}{n}\sum_{i=1}^{n}\left(\frac{\left(Y_{i}-f_{1}(W_{i})\right)D_{i}}{p(W_{i})}-\frac{\left(Y_{i}-f_{0}(W_{i})\right)(1-D_{i})}{1-p(W_{i})}\right).
\end{align*}
$$
This is very clever construction because it has two equivalent forms: ECATE+(roughly) mean zero error and IPW+(roughly) mean zero error:

- In the current form, the first part is exactly $\hat{\tau}^{*}_{ECATE}$. The second part, instead of an estimator of ATE, is actually a weighted sum of two mean zero terms. For instance, by unconfoundedness,
  $$
  \begin{align*}
  &\mathbb{E}\left[\left(Y_{i}-f_{1}(W_{i})\right)D_{i}\right]\\
  =&\mathbb{E}\left[\mathbb{E}\left(\left(Y_{i}-f_{1}(W_{i})\right)D_{i}|W_{i}\right)\right]\\
  =&\mathbb{E}\left[\Pr(D_{i}=1|W_{i})\mathbb{E}\left(\left(Y_{i}-f_{1}(W_{i})\right)D_{i}|D_{i}=1,W_{i}\right)\right]\\
  =&\mathbb{E}\left[\Pr(D_{i}=1|W_{i})\mathbb{E}\left(\left(Y_{i}-f_{1}(W_{i})\right)|D_{i}=1,W_{i}\right)\right]\\
  =&\mathbb{E}\left[\Pr(D_{i}=1|W_{i})\left(f_{1}(W_{i})-f_{1}(W_{i})\right)\right]\\
  =&0.
  \end{align*}
  $$

  - This suggests that, later on, when we need to estimate the **propensity score**, its estimation error would not much affect the estimator of ATE because the estimated propensity score inverse is multiplied by a mean zero term, leading to ignorable asymptotic bias. So **the propensity score NEEDS NOT to be consistently estimated as long as the conditional expectations $f_{d}$ are consistently estimated**.

- Now here is the real magic: Like a transformer, $\hat{\tau}^{*}_{AIPW}$ has a second form. By simple algebra, we can transform $\hat{\tau}^{*}_{AIPW}$ into the following:
  $$
  \begin{align*}
  \hat{\tau}^{*}_{AIPW}=&\frac{1}{n}\sum_{i=1}^{n}\left(\frac{Y_{i}D_{i}}{p(W_{i})}-\frac{Y_{i}(1-D_{i})}{1-p(W_{i})}\right)\\
  &+\frac{1}{n}\sum_{i=1}^{n}\left(f_{1}(W_{i})\left(1-\frac{D_{i}}{p(W_{i})}\right)-f_{0}(W_{i})\left(1-\frac{1-D_{i}}{1-p(W_{i})}\right)\right).
  \end{align*}
  $$

  - The first term is exactly $\hat{\tau}_{IPW}^{*}$! The second term, instead of an estimator of ATE, is again a weighted sum of two mean zero terms. For instance,

  $$
  \mathbb{E}\left[1-\frac{D_{i}}{p(W_{i})}\right]=1-\mathbb{E}\left[\frac{1}{p(W_{i})}\mathbb{E}\left(D_{i}|W_{i}\right)\right]=0.
  $$

  - This suggests that, later on, when we need to estimate the **conditional expectations**, their estimation errors would not much affect the estimator of ATE because the estimated conditional expectations are multiplied by mean zero terms, leading to ignorable asymptotic bias. So **the conditional expectations NEED NOT to be consistently estimated as long as the propensity score function is consistently estimated**.

### 3.2 Weak Double Robustness

Let $\hat{f}_{d}$ and $\hat{p}$ be estimators of $f_{d}$ and $p$. Now we can define a feasible estimator of ATE as follows:
$$
\begin{align*}
\hat{\tau}_{AIPW}\equiv&\frac{1}{n}\sum_{i=1}^{n}\left(\hat{f}_{1}(W_{i})-\hat{f}_{0}(W_{i})\right)\\
&+\frac{1}{n}\sum_{i=1}^{n}\left(\frac{\left[Y_{i}-\hat{f}_{1}(W_{i})\right]D_{i}}{\hat{p}(W_{i})}-\frac{\left[Y_{i}-\hat{f}_{0}(W_{i})\right](1-D_{i})}{1-\hat{p}(W_{i})}\right).
\end{align*}
$$
By the analysis at the end of Section 3.1, we know that $\hat{\tau}_{AIPW}$ is **weakly doubly robust**: $\hat{\tau}_{AIPW}\overset{p}{\to}\tau\equiv ATE$ as long as one of $(\hat{f}_{1},\hat{f}_{0})$ and $\hat{p}$ is consistent!

This is a very interesting result. It says that to consistently estimate ATE under unconfoundedness, you only need to correctly estimate either the CATE function or the propensity score, and, most remarkably, you don't need to know which one is wrong. This is in sharp contrast to classical statistical/econometric methods where any mistake could result in a wrong final estimate.

- If you heard about GMM, you know we can in fact estimate $\tau$ by GMM since we have two moment equations ((3) and (7)) but only one parameter. However, GMM would be inconsistent if any of those two moment estimators is inconsistent.

### 3.3 Strong Double Robustness and Double Machine Learning

Weak double robustness only tells us consistency. We have a few concerns about the other aspects about $\hat{\tau}_{AIPW}$.

- Unlike the infeasible AIPW $\hat{\tau}^{*}_{AIPW}$ where the estimator is simply constructed by some averages of i.i.d. random variables, the feasible version $\hat{\tau}_{AIPW}$ is constructed by the average of some estimated functions. Chapter 6 tells us that these estimators, no matter whether following the classical nonparametrics or the modern machine learning methods,  have a variance whose order is larger than $1/n$, resulting in a large standard error and a **slow rate of convergence** (i.e., although the estimator is consistent, it converges in probability to the true ATE slowly).
- In data sciences, we usually have the belief of "garbage in, garbage out". Now that the input has a slow rate of convergence, one may expect to see that the final estimator also has a slow rate of convergence. Meanwhile, the variance of the final estimator might also be large because the variances of the nonparametrically or ML estimated $f_{d}$ and $p$ are large.

However, in this section, we introduce a modified estimator based on our $\hat{\tau}_{AIPW}$ which achieves **strong double robustness**: The estimator's variance **has order 1/n** as long as the order of the variances of $\hat{f}_{d}$ and $\hat{p}$ are **not too large**. In other words, the estimator's rate of convergence achieves **the fastest speed** (i.e., the speed of $\hat{\tau}^{*}_{AIPW}$) as long as the convergence speed of $\hat{f}_{d}$ and $\hat{p}$ is **not too slow**.  In particular, **causal tree** and **causal forest** estimators of $f_{d}$ and $p$ satisfy the requirements.

We now introduce this modified estimator. It is an example of **double machine learning** (DML), which uses the idea of **honest splitting** in causal trees. Recall that in the infeasible estimator $\hat{\tau}_{AIPW}^{*}$, we use data (more specifically the $W_{i}$s) to estimate ATE. However, in the feasible version $\hat{\tau}_{AIPW}$, the data are used twice for different purposes; they are first used to estimate the conditional expectations and propensity score, and then employed to construct the final estimate. Recall that we had a similar issue in the causal tree/forest estimator. So we can do the following honest splitting and **cross-fitting**:

- Split your data set randomly into two halves $\mathcal{I}_{1}$ and $\mathcal{I}_{2}$.

- Use data in $\mathcal{I}_{2}$ to estimate $f_{d}$ and $p$, by, for example, causal forest. Call them $\hat{f}_{d}^{\mathcal{I}_{2}}$ and $\hat{p}^{\mathcal{I}_{2}}$.

- Plug in $W_{i}$ in $\mathcal{I}_{1}$ into them and construct the following estimator. 
  $$
  \begin{align*}
  \hat{\tau}^{\mathcal{I}_{1}}\equiv&\frac{1}{|\mathcal{I}_{1}|}\sum_{i\in\mathcal{I}_{1}}\left(\hat{f}_{1}^{\mathcal{I}_{2}}(W_{i})-\hat{f}_{0}^{\mathcal{I}_{2}}(W_{i})\right)\\
  &+\frac{1}{|\mathcal{I}_{1}|}\sum_{i\in\mathcal{I}_{1}}\left(\frac{\left[Y_{i}-\hat{f}_{1}^{\mathcal{I}_{2}}(W_{i})\right]D_{i}}{\hat{p}^{\mathcal{I}_{2}}(W_{i})}-\frac{\left[Y_{i}-\hat{f}_{0}^{\mathcal{I}_{2}}(W_{i})\right](1-D_{i})}{1-\hat{p}^{\mathcal{I}_{2}}(W_{i})}\right).
  \end{align*}
  $$

- Now, swap the roles of $\mathcal{I}_{1}$ and $\mathcal{I}_{2}$ and construct $\hat{\tau}^{\mathcal{I}_{2}}$. 

  - This step is called **cross-fitting**. It improves the stability of the final estimator because every data point is eventually used to do to jobs, estimating the nonparametric functions and estimating the ATE.

- Construct our double machine learning estimator of ATE as
  $$
  \hat{\tau}^{DML}_{AIPW}\equiv \frac{|\mathcal{I}_{1}|}{n}\hat{\tau}^{\mathcal{I}_{1}}+\frac{|\mathcal{I}_{2}|}{n}\hat{\tau}^{\mathcal{I}_{2}}.
  $$

It turns out $\hat{\tau}^{DML}_{AIPW}$  has a lot of nice properties under two assumptions: unconfoundedness and **strong overlap**. Strong overlap assumes that $p(W_{i})$ is not too close to 1 and 0 for all realizations of $W_{i}$. It means for every value of $W$, the proportion of people selecting either $D=1$ or $D=0$ is not too small. For instance, if $W_{i}$ means gender and $D_{i}$ means whether to get a master's degree, then strong overlap says that for each gender, there are not too few people who have a master's degree and not too few people who do not have one. 

Under unconfoundedness and strong overlap:

- $\hat{\tau}^{DML}_{AIPW}$ is both weakly and strongly doubly robust. In particular, let the order of variance of $\hat{f}_{d}$ be $1/n^{a}$ and that of $p$ be $1/n^{b}$. Then the variance of $\hat{\tau}^{DML}_{AIPW}$ has order $1/n$ as long as $a+b\geq 1$, 

- Moreover, there is a positive number $V_{AIPW}$ such that
  $$
  \sqrt{n}\left(\hat{\tau}^{DML}_{AIPW}-\tau\right)\overset{d}{\to}\mathcal{N}(0,V_{AIPW}).
  $$

- More impressively, among **a very large class of estimators** of ATE, $\hat{\tau}^{DML}_{AIPW}$ has **the smallest variance**.

- To make things even better, $V_{AIPW}$ can be consistently estimated, so the standard error $se$ can be easily obtained by $\sqrt{\hat{V}_{AIPW}/n}$. Then, $t$-test, p-value and confidence interval can be constructed in the standard way.

- To sum up, $\hat{\tau}^{DML}_{AIPW}$ is the fastest, most efficient, and most robust estimator of ATE under unconfoundedness. In short, it is the best. Always do it as long as you believe in you data unconfoundedness holds but you don't have complete randomization.

## 4. Implementation

Because of the strong double robustness, you can choose most nonparametric/machine learning methods to estimate $f_{d}$ and $p$. Plugging them into the formula of $\hat{\tau}^{DML}_{AIPW}$ yields the desired estimator. For convenience, we recommend to estimate $f_{d}$ and $p$ by causal forest because one simple command **average_treatment_effect** in the grf package produces the estimated ATE and the standard error.

```R
library(grf)
library(DiagrammeR)
### Generate data. 
set.seed(5280)
n <- 2000      
p <- 10
W <- matrix(rnorm(n * p), n, p)
D <- rbinom(n, 1, 0.4 + 0.2 * (W[, 1] > 0))
Y <- pmax(W[, 1], 0) * D + W[, 2] + pmin(W[, 3], 0) + rnorm(n)
data=data.frame(Y,D,W)

# Plant a forest
tau.forest <- causal_forest(W, Y, D, num.trees = 4000)

# Estimate the ATE
average_treatment_effect(tau.forest,target.sample="all",method="AIPW")
```











 

