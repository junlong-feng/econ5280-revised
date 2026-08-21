# ECON5280 Chapter 8 Sharp Regression Discontinuity Design

<font size="5">Junlong Feng</font>

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/junlong-feng/econ5280/main?filepath=Lecture8_RD.ipynb)

## Outline

* Motivation: In some very special cases, $D$ is a deterministic function of $W$.
* Identification of ATE or CATE: ATE is impossible whereas CATE can only be identified at ONE value of $W$.
* Locally Linear Estimation: How to, pros, and cons.
* Applications and Implementation.

## 1. A Very Special Link Between $D$ and $W$

In the past two chapters, we've been working under the unconfoundedness assumption:
$$
D_{i}\perp (Y_{i}(1),Y_{i}(0))|W_{i}.
$$
This assumption says that given $W$, the random variable $D$ is independent of the potential outcomes. There is one very special case of unconfoundedness, that is, $D_{i}$ is a **deterministic function** of $W_{i}$. In that case, conditional on $W_{i}$, $D_{i}$ is **nonrandom**; independence trivially holds, but it loses most of its identification power. To see it, recall that for any value of $W_{i}$, denoted by $w$, we have the following identication equation **provided that the two conditional expectations on the right hand side is well defined**:
$$
CATE(w)\overset{unconfoundedness}{=}\mathbb{E}\left(Y_{i}|D_{i}=1,W_{i}=w\right)-\mathbb{E}\left(Y_{i}|D_{i}=0,W_{i}=w\right).
$$
Let's try to figure out what the problem is if $D$ is a deterministic function of $W$ using the following example.

- Suppose $W_{i}$ is a scalar and $D_{i}=1$ if and only if $W_{i}\geq w_{0}$.
- For instance, suppose there is a policy saying that pregnant women whose (family) income $(W_{i})$ is below $w_{0}$ can receive free maternal care $(D_{i}=0)$. 
- If $w>w_{0}$, then $\mathbb{E}(Y_{i}|D_{i}=0,W_{i}=w)$ is not well defined: It is impossible for anyone whose income is above the threshold but receive free maternal care.
- Similarly, if $w<w_{0}$, then $\mathbb{E}(Y_{i}|D_{i}=1,W_{i}=w)$ is not well defined.
- Hence, identification equation (2) fails in general.

In this case, it is impossible to identify the entire CATE function. Notice that there's no need to specify $D$ now: $D_{i}=1$ if and only if $W_{i}\geq w_{0}$.  Therefore, by the definition of CATE: 
$$
\begin{equation}
\begin{aligned}
CATE(w)\equiv&\mathbb{E}(Y_{i}(1)|W_{i}=w)-\mathbb{E}(Y_{i}(0)|W_{i}=w)\\
=&\begin{cases}
\mathbb{E}(Y_{i}|W_{i}=w)-\mathbb{E}(Y_{i}(0)|W_{i}=w),&\text{if }w\geq w_{0};\\
\mathbb{E}(Y_{i}(1)|W_{i}=w)-\mathbb{E}(Y_{i}|W_{i}=w),&\text{if }w<w_{0}.
\end{cases}
\end{aligned}
\end{equation}
$$
So for either $w>w_{0}$ or $w<w_{0}$, there is one conditional expectations of the potential outcome that cannot be identified. This prevents us from identifying $CATE(w)$ for almost every $w$.

## 2. Identification of Regression Discontinuity Design

However, there is a very special point: $W_{i}=w_{0}$. 

- Suppose both $\mathbb{E}(Y_{i}(1)|W_{i}=w)$ and $\mathbb{E}(Y_{i}(0)|W_{i}=w)$ are continuous functions in $w$ at $w_{0}$. 

- Continuity and equation (3) imply that
  $$
  \begin{align*}
  \mathbb{E}(Y_{i}(1)|W_{i}=w)=\lim_{w\downarrow w_{0}}\mathbb{E}(Y_{i}|W_{i}=w),\\
  \mathbb{E}(Y_{i}(0)|W_{i}=w)=\lim_{w\uparrow w_{0}}\mathbb{E}(Y_{i}|W_{i}=w).
  \end{align*}
  $$
  See the class notes for a graphical illustration.

- This implies that $CATE(w_{0})$ is identified:
  $$
  CATE(w_{0})=\lim_{w\downarrow w_{0}}\mathbb{E}(Y_{i}|W_{i}=w)-\lim_{w\uparrow w_{0}}\mathbb{E}(Y_{i}|W_{i}=w).
  $$
  
- Variable $W_{i}$ is often called the **running variable**.

Identification here is **NOT DUE TO** unconfoundedness.

- **Two key components** for identification: i) The deterministic relationship between $D$ and $W$ which yields equation (3), and ii) continuity in $\mathbb{E}(Y_{i}(d)|W_{i}=w)$ at $w_{0}$ for both $d=1$ and $0$.
- The treatment effect is equal to the jump from the left to the right limit of $\mathbb{E}(Y_{i}|W_{i}=w)$ at $w_{0}$; If there is no jump, then the treatment effect is 0. 
- Identification by this type of discontinuity of $\mathbb{E}(Y_{i}|W_{i}=w)$ at $w_{0}$ is called **regression discontinuity design** (RDD).

Some economic intuition of the source of identification:

- Consider the maternal care policy as an example. If you want to study the causal effect of this policy on children's birthweight (birthweight is an important health indicator which has long run causal effect on a child's future cognitive development), then essentially you want to compare the birthweight of the children whose mothers got the free maternal care with that of children whose mothers did not. However, by the policy design, these two groups of mothers had very different financial situations when they were pregnant: mothers who got the free care were systematically poorer than those who did not. Since mothers' or family financial status may also causally impact children's birthweight, the proposed comparison is biased by it.
- Now instead, if we only focus on the causal effect of this policy on a very special group of children: children whose mothers' (family) income was at (or very close to) the policy cutoff line $w_{0}$. For example, suppose $w_{0}=\$1000$. Those who earn $\$1001$ a month may not be systematically different from those who earn $\$999$ a month. This two-dollar difference may simply be driven by random stuff. However, in terms of the treatment assignment, women at these two salary levels are assigned to two totally different groups. The difference of these two groups' children's birthweight can be convincingly credited to the free maternal care. Causal effect for these people is identified.

No matter from the mathematical derivation or from the intuition, the identified CATE is **just one value** $CATE(w_{0})$ rather than **a whole function** $CATE(\cdot)$ as in Chapter 6. Note that this single value **may or may not** be interesting.

- Suppose the policy really wants to help those extremely low income pregnant women, i.e., $W_{i}<<w_{0}$, say monthly income equal to $\$500$. 
- Then the identified causal effect **fails to answer** the question whether the policy is useful to them because we only know the effectiveness of the policy for those with $W_{i}=\$1000$. 

Therefore, you should be very careful when you're using and interpreting this estimate. In general, you cannot generalize your findings too much. They may be not extendable to other groups of individuals.

## 3. Locally Linear Estimation

According to the identification equation (4), to estimate $CATE(w_{0})$, we need to estimate $E(Y_{i}|W_{i}=w)$ for $w<w_{0}$ and $w>w_{0}$, and take the left and right limits of the estimators, respectively. So far, we have learned tree methods to estimate conditional expectations. In practice however, economists like to estimate them by **locally linear regression**, a standard local nonparametric method which is very easy to handle.

Recall from calculus that for sufficiently smooth functions, a small chunk of it, no matter how curvy, can always be well-approximated by straight lines.  So this is our plan:

- Pick a small neighborhood around $w_{0}$: $[w_{0}-h,w_{0}+h]$. 

- Now we are going to fit two straightlines, one to the left of $w_{0}$ and one to the right. They can of course have different slopes and intercepts. So let them be:
  $$
  \begin{align*}
  Y_{i}=&\beta_{0}^{+}+\beta_{1}^{+}(W_{i}-w_{0})+\varepsilon^{+}_{i},\ \ \forall W_{i}\geq w_{0},\\
  Y_{i}=&\beta_{0}^{-}+\beta_{1}^{-}(W_{i}-w_{0})+\varepsilon^{-}_{i},\ \ \forall W_{i}<w_{0}.
  \end{align*}
  $$
  Then we have
  $$
  CATE(w_{0})\approx \beta_{0}^{+}-\beta_{0}^{-}.
  $$

- However, we don't have to run two regressions to separately obtain $\beta_{0}^{+}$ and $\beta_{0}^{-}$. Let $D_{i}=1(W_{i}\geq w_{0})$ and write
  $$
  \begin{align*}
  Y_{i}=&\gamma_{0}+\gamma_{1}D_{i}+\gamma_{2}(W_{i}-w_{0})+\gamma_{3}D_{i}(W_{i}-w_{0})+\varepsilon_{i}.
  \end{align*}
  $$
  
- One can verify that

  - $\gamma_{0}=\beta_{0}^{-}$.
  - $\gamma_{2}=\beta_{1}^{-}$.
  - $\gamma_{0}+\gamma_{1}=\beta_{0}^{+}$.
  - $\gamma_{2}+\gamma_{3}=\beta_{1}^{+}$.

  So, $\gamma_{1}$ is the $CATE(w_{0})$.

- Pick a weighting function $k$ such that higher weights will be given to data points whose $W_{i}$ are closer to $w_{0}$. Function $k$ is positive, symmetric around 0, and attains its unique maximum at 0. So, the weights are equal to $k\left(\frac{W_{i}-w_{0}}{h}\right)$.

- Let $X_{i}\equiv(1,D_{i},W_{i}-w_{0},D_{i}(W_{i}-w_{0}))'$. Let $\gamma\equiv (\gamma_{0},\gamma_{1},\gamma_{2},\gamma_{3})'$. Estimate $\gamma$ by the weighted method of moments (or equivalently, weighted OLS):
  $$
  \begin{align*}
  &\frac{1}{n}\sum_{i=1}^{n}k\left(\frac{W_{i}-w_{0}}{h}\right)X_{i}\left(Y_{i}-X_{i}'\hat{\gamma}\right)=0,\\
  \implies& \hat{\gamma}=\left(\frac{1}{n}\sum_{i}k\left(\frac{W_{i}-w_{0}}{h}\right)X_{i}X_{i}'\right)^{-1}\left(\frac{1}{n}\sum_{i}k\left(\frac{W_{i}-w_{0}}{h}\right)X_{i}Y_{i}\right).
  \end{align*}
  $$
  
- 

This locally linear estimator of $\gamma_{1}$, i.e., the second entry in $\hat{\gamma}$, is consistent of $CATE(w_{0})$, and asymptotically normal with easy-to-calculate standard errors if $h\to 0$ with $n\to\infty$. So, testing, p-value, and confidence interval can be obtained in the standard way.

- But, of course, due to its nonparametric nature, its variance is relatively large and it converges to the true CATE in probability slowly since we only used a subsample.

### 3.1 Bandwidth Selection

In theory, we need $h\to 0$ at a proper rate. However, in practice, we need to pick a concrete number. An often-adopted bandwidth selection is called **IK bandwidth** by Imbens and Kalyanaraman (2012, RES). It is computed by a complicated algorithm, but standard R packages will compute it for you. Packages such as "rdd" also reports estimates based on half and twice of the IK bandwidth to show robustness.

### 3.2 Including Covariates

So far, our $W_{i}$ is a scalar which determines $D_{i}$. There is **no need to include other covariates for identification** if we do not care about the CATE in terms of them. One motivation to include covariates is to reduce the standard error like in Chapter 5. If that is the goal, applied researchers usually add the covariates linearly into the locally linear regression (Calonico, Cattaneo, Farrell, and Titianik, 2012, ReStat). More work needed.

### 3.3 Alternative Estimators

Recall that we can also use **global methods** to approximate a function. Indeed, an often-adopted estimator besides the locally linear estimator is global polynomial approximation, which we explained in Chapter 6. The empirical example we'll see in the next section uses it. 

## 4. Applications and Implementation

Regression discontinuity design has a wide application in economics, political science, and other social sciences. In particular, it’s a workhorse for electoral research.

- David Lee (2008, Journal of Econometrics; the figure below is taken from it) studies the effect of incumbency on the chance of winning an election.
- Hypothesis: Incumbency may have impact on voter’s voting decision, so may increase the chance of winning the election.
- Data: District level data for U.S. house election. Democratic vs. Republican. 
- Naïve idea: regress the vote share for democratic on whether the vote share was over 50% for democratic in previous election.
- Endogeneity: Winning the previous election may be because democrats are good, which may lead to another win.
- RDD: look at districts where the vote share for Democratic in the previous election was just around 50%. So winning or losing was almost random. Not much to do with whether Democrats were really doing good.

<img src="/Users/junlong/Downloads/Ch8RD_Lee.png" style="zoom:50%;" />

- Running variable: Democratic vote share (recentered around 0.5) in election $t$.
- Outcome variable: Democratic vote share in election $t+1$.
- You can see that these two vote shares are positively correlated no matter whether Democrats won last time.
- But winning in the last time does cause a jump in the vote share this time.
- This jump is the incumbency effect.
- The dots are averages using data lying close to each other; not plotting all the original data just because it’s too messy.
- The curves are fitted using polynomials.

Other applications:

- Wage increase on performance of mayors (Ferraz and Finan 2011; Gagliarducci and Nannicini 2013)

- Colonial institutions on development outcomes (Dell 2009)

- Financial aid offers on college enrollment (Van der Klaauw 2002)

- Access to Angel funding on growth of start-ups (Kerr, Lerner, and Schoar 2010).

### 4.1 Implementation

There are more than one R packages to do RDD. We chose a package called "rdd" which runs locally linear regression for estimation. The following example is taken from its documentation.

```R
library(rdd)

## DGP
x<-runif(1000,-1,1)
y<-3+2*x+10*(x>=0)+rnorm(1000)

## Estimation; the default bandwidth is IK and kernel is triangular. You can specify them by "bw=" and "kernel="
## If you want to add covariates (cov), then write RDestimate(y~x|cov). 
model1=RDestimate(y~x) 
model2=RDestimate(y~x,bw=c(0.1,0.2,0.3),kernel="gaussian")
summary(model1)
summary(model2)
plot(model1)
```

Recall that Lee's results are based on polynomial regression. We can estimate it by locally linear regression and see whether the results are robust:

```R
library(rdd)
data=read.csv("house.csv")
y=data$y
x=data$x
model=RDestimate(y~x)
summary(model)
```

