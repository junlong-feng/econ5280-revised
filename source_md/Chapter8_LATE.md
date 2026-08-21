# ECON5280 Chapter 8 Instrumental Variable

<font size="5">Junlong Feng</font>

Outline

* Motivation: What shall we do when there are compliance issues?
* Noncompliance and LATE: Average treament effect for a subpopulation.
* Identification of LATE: When is LATE identified (in a linear model)?
* Estimation: Another MM estimator.
* Applications: How to find an instrument in practice.

## 1. Noncompliance and LATE (Binary Treament and Binary IV)

Sometimes we can guarantee the treatment assignment is fully random, but not everyone comply.

- Suppose the government freely assign HPV vaccine to citizens by a lottery. If Alice wins the lottery, she can get vaccinated for free.
- This lottery is fully random and independent of everything else.
- Suppose I'm a health economist and want to study the effect of receiving HPV vaccine on newborns' birthweight $Y$.
- Can I use this lottery as $Z$ and estimate $Y=\beta_{0}+\beta_{1}Z+\varepsilon$ by OLS?
- No because:
  - i) Winning the lottery DOES NOT imply that the individual must go and get the vaccine; people who don't want to get vaccinated won't go anyways.
  - ii) Not winning the lottery DOES NOT imply the individual cannot get the vaccine; people who want to get vaccianted can pay for it.
  - So $\beta_{1}$ is the ATE of winning the lottery on birthweight, not the ATE of getting HPV vaccine on birthweight. 

- This ATE is called intent-to-treat (ITT).

So can we get ATE? No because whether people get vaccinated ($D$) is determined by a lot of factors, which are in the error $\varepsilon$. However, we can get the average effect of a sub-population.

- Suppose Alice is on the margin: She's quite indifferent whether to get the vaccine or not.
- This means that the marginal cost of getting vaccination (money, time, etc) is equal to her marginal benefit.
- Now suppose Alice wins the lottery, the marginal cost becomes lower, and she decides to get vaccinated.
- This decision $D=1$ is made solely based on winning the lottery, independent of everything else.
- So for Alice, $D$ is as good as random.
- Now suppose Bob doesn't want to get vaccinated anyways. He has already made up his mind based on his $\varepsilon$. No matter whether his $Z=1$ or not, he won't do it. So his decision $D=1$ is correlated with $\varepsilon$, not randomly assigned.

The above example shows that it's useful to divide the population into subgroups more carefully. Define $D(z)$ as the potential treatment. $D=D(z)$ if and only if $Z=z$. 

The above example shows that we can divide the population into several groups 

| Always Taker  |  Never Taker  |  Complier   |   Defier    |
| :-----------: | :-----------: | :---------: | :---------: |
| $D(1)=D(0)=1$ | $D(1)=D(0)=0$ | $D(1)>D(0)$ | $D(1)<D(0)$ |

- We never know which group a given individual $i$ belongs to because we can never observe her both potential treatments. Same logic as the *fundamental problem of causal inference*.
- Compliers and defiers are relative; you can switch their definitions according to applications. In the HPV example, it's more natural to define complier in the above way.

From the vaccination example, we learned that the randomly assigned $Z$, which is called an **instrumental variable, or instrument, or IV**, can represent $D$ for compliers. So we may conjecture that the ATE for compliers might be identified. Let's first define the ATE for this subgroup:

**Definition (LATE)**. The local average treatment effect (LATE) is defined as 
$$
LATE\equiv \mathbb{E}\left(Y_{i}(1)-Y_{i}(0)|D_{i}(1)>D_{i}(0)\right).
$$
where the conditioning part means that $i\in Compliers$.

* Without the conditioning part, the RHS is ATE.
* If everyone in the population is a complier, then $D_{i}(1)>D_{i}(0)$ holds with probability 1 so $LATE=ATE$.

## 2. Identification of LATE (Binary Treament and Binary IV)

We make the following assumptions:

- $Z\perp (Y(1),Y(0),D(1),D(0))$. 
  - This assumption says that the instrument is completely randomized. 
  - Also called **exclusion**.
  - We can extend it to conditional randomized as in Chapter 6. Will do later.
- $\mathbb{E}(D|Z=1)\neq \mathbb{E}(D|Z=0)$. 
  - This assumption says not all individuals are always takers or never takes. Because otherwise, $D=1$ or 0 with probability 1 so $\mathbb{E}(D|Z=1)=\mathbb{E}(D|Z=1)$.
  - Also called **relevance**.
- $D(1)\geq D(0)$ with probability 1.
  - This means there are no defiers in the population.
  - Also called **monotonicity**.

**Theorem (LATE, Imbens and Angrist, 1994)**. Under the above assumptions, LATE is identified as
$$
LATE=\frac{\mathbb{E}(Y|Z=1)-\mathbb{E}(Y|Z=0)}{\mathbb{E}(D|Z=1)-\mathbb{E}(D|Z=0)}.
$$

- This result is called **identification** because the unknown parameter of interest LATE is uniquely linked to some population quantity that is directly estimable.

**Proof**. Recall that $D=D(0)$ if $Z=0$ and $D=D(1)$ if $Z=1$, where $D(0),D(1)$ are potential treatment. Equivalently, $D=ZD(1)+(1-Z)D(0)$. Meanwhile, $Y=DY(1)+(1-D)Y(0)$. Therefore,
$$
\begin{align*}
\mathbb{E}(Y|Z=1)-\mathbb{E}(Y|Z=0)=&\mathbb{E}\left[D(1)Y(1)+(1-D(1))\cdot Y(0)|Z=1\right]\\
&-\mathbb{E}\left[D(0)Y(1)+(1-D(0))\cdot Y(0)|Z=0\right]\\
\text{By independence:}=&\mathbb{E}\left[D(1)Y(1)+(1-D(1))\cdot Y(0)\right]\\
&-\mathbb{E}\left[D(0)Y(1)+(1-D(0))\cdot Y(0)\right]\\
=&\mathbb{E}\left[\left(D(1)-D(0)\right)\times \left(Y(1)- Y(0)\right)\right]\\
\text{By monotonicity: }=&\mathbb{E}\left[Y(1)- Y(0)|D(1)>D(0)\right]\times \Pr(D(1)-D(0)=1)\\
\text{By independence: }=&\mathbb{E}\left[Y(1)- Y(0)|D(1)>D(0)\right]\times \left[\mathbb{E}(D|Z=1)-\mathbb{E}(D|Z=0)\right].
\end{align*}
$$
Done by dividing $\mathbb{E}(D|Z=1)-\mathbb{E}(D|Z=0)$ on both side because it's nonzero by the relevance condition. 

### 2.1 Linear Representation

The identification result for LATE looks formidable: a lot of expectations. It is possible to have a linear representation for it, just as ATE under complete randomization.

- $\gamma\equiv \mathbb{E}(D|Z=1)-\mathbb{E}(D|Z=0)$. By identity
  $$
  D=\mathbb{E}(D|Z=0)+[\mathbb{E}(D|Z=1)-\mathbb{E}(D|Z=0)] \times Z+(D-\mathbb{E}(D|Z)).
  $$
  We have
  $$
  D=\gamma_{0}+\gamma Z+\nu,
  $$
  where $cov(Z,\nu)=0$ is ALWAYS true without any assumption. (Can you prove it?)

  - By assumption $\gamma\neq 0$.

- $\delta\equiv \mathbb{E}(Y|Z=1)-\mathbb{E}(Y|Z=0)$. Similarly, we can write
  $$
  Y=\delta_{0}+\delta Z+\mu,
  $$
  where $\delta_{0}=\mathbb{E}(Y|Z=0)$ and $\mu=Y-\mathbb{E}(Y|Z)$ so $cov(Z,\mu)=0$ without any assumptions as well.

- Combining them, let $\beta=LATE=\delta/\gamma$ (this last equality is **by the LATE theorem**), we have
  $$
  \begin{align*}
  Y=&\delta_{0}+\delta\frac{D-\gamma_{0}-\nu}{\gamma}+\mu\\
  =&\left(\delta_{0}-\delta\frac{\gamma_{0}}{\gamma}\right)+\frac{\delta}{\gamma}D+\left(\mu-\frac{\delta}{\gamma}\nu\right)\\
  \eqqcolon&\beta_{0}+\beta_{1} D+\varepsilon.
  \end{align*}
  $$
  
- **Important**. $cov(D,\varepsilon)\neq 0$ because $D$ is related to $\nu$, but $cov(\varepsilon,Z)=0$.

**Theorem**. When $D$ and $Z$ are binary and under the assumptions for the LATE theorem, there exists a unique linear model $Y_{i}=\beta_{0}+\beta_{1} D_{i}+\varepsilon_{i}$ such that i) $\beta_{1}=LATE$ and ii) $cov(Z_{i},\varepsilon_{i})=0$ and $\mathbb{E}(\varepsilon_{i})=0$.

### 2.2 Conditional LATE (CLATE)

Similar to ATE and CATE, we can relax the assumptions for the LATE theorem to hold conditionally. Suppose we have a vector of control variables $W$ such that for $Y=g(D,W,U)$ and $D=h(D,W,V)$, 

- $Z\perp (U,V)|W$. Or equivalently, $Z\perp (Y(1),Y(0),D(1),D(0))|W$. 
  - This assumption says the instrument is conditional randomized. 
- $\mathbb{E}(D|Z=1,W)\neq \mathbb{E}(D|Z=0,W)$. 
- $D(1)\geq D(0)$ with probability 1 conditional on $W$.

Then defining the conditional LATE (CLATE) as 
$$
CLATE(W)\equiv \mathbb{E}\left(Y_{i}(1)-Y_{i}(0)|D_{i}(1)>D_{i}(0),W\right).
$$
We can show that 
$$
CLATE(W)=\frac{\mathbb{E}(Y|Z=1,W)-\mathbb{E}(Y|Z=0,W)}{\mathbb{E}(D|Z=1,W)-\mathbb{E}(D|Z=0,W)}.
$$
**Important**. CLATE no longer always has a representation in linear model, just like CATE. One can estimate it by instrument forest using the grf R package.

## 3 Instrumental Variable Regression

By the linear representation, estimating LATE boils down to estimating $\beta_{1}$ in
$$
Y_{i}=\beta_{0}+\beta_{1} D_{i}+\varepsilon_{i},\ \ \mathbb{E}(\varepsilon_{i}Z_{i})=0\ \&\ \mathbb{E}(\varepsilon_{i})=0.
$$
We can estimate the coefficients by an MM estimator:

- Moment equations:
  $$
  \begin{align*}
  \mathbb{E}(Y_{i}-\beta_{0}-\beta_{1}D_{i})&=0,\\
  \mathbb{E}[(Y_{i}-\beta_{0}-\beta_{1}D_{i})Z_{i}]&=0.
  \end{align*}
  $$

- The first equation says $\beta_{0}=\mathbb{E}(Y_{i})-\beta_{1}\mathbb{E}(D_{i})\equiv \mu_{Y}-\beta_{1}\mu_{D}$.

- Substitute it into the second:
  $$
  \mathbb{E}[(Y_{i}-\mu_{Y})(Z_{i}-\mu_{Z})]=\beta_{1}\mathbb{E}[(D_{i}-\mu_{D})(Z_{i}-\mu_{Z})].
  $$

- Therefore, $LATE=\beta_{1}=cov(Y_{i},Z_{i})/cov(D_{i},Z_{i})$, **provided that $cov(D_{i},Z_{i})\neq 0$**.

  - This condition is implied by the relevance condition because
    $$
    \begin{align*}
    cov(D_{i},Z_{i})=&\mathbb{E}(D_{i}Z_{i})-\mu_{D}\mu_{Z}\\
    =&\mathbb{E}(D_{i}|Z_{i}=1)\Pr(Z_{i}=1)-\mu_{D}\mu_{Z}\\
    =&\mu_{Z}(\mathbb{E}(D_{i}|Z_{i}=0)-\mathbb{E}(D_{i})),
    \end{align*}
    $$
    which is nonzero if and only if $\mathbb{E}(D_{i}|Z_{i}=0)\neq \mathbb{E}(D_{i}|Z_{i}=1)$.

- We can estimate $\beta_{1}$ by replacing all the expectations with sample averages:
  $$
  \hat{\beta}^{IV}\equiv\frac{\sum_{i=1}^{n}(Y_{i}-\bar{Y})(Z_{i}-\bar{Z})}{\sum_{i=1}^{n}(D_{i}-\bar{D})(Z_{i}-\bar{Z})}.
  $$

$\hat{\beta}^{IV}$ is consistent and asymptotically normal. Inference can be done by standard methods.

### 3.1 Relevance: Weak Instrument and First Stage F-value

When $Z$ and $D$ are only weakly correlated, $Z$ is a **weak instruments**. Weak IVs cause a lot of problems, e.g. bias, poor normal approximation of the distribution, etc. See Q.4 in Problem Set 3 for an example.

There are formal test to run, but when you only have one endogenous variable, you can regress the endogenous variable on $Z$, and look at the F-value of that regression. 

- If the F-value is greater than 10, it suggests the IV are not weak.

Weak IV is still an active research area.

## 4. Applications

Instrumental variable regression is the workhorse in economics. It revolutionizes the entire discipline, making causal inference possible for social sciences in many important areas where it's impossible to run experiment.

- Two decades ago you get a PhD degree and a job in a nice econ department if you find a smart IV.
- IV is still one of the most popular method in applied econ, and on top of that, some other equally popular method, e.g., fuzzy regression discontinuity design (RDD), are just special IVs in particular applications.

Implementation is very easy and thus omitted. You can use command ivreg from the R package "ivreg" to do everything. A more difficult problem in practice is where to find an IV.

### 4.1 IVs from RCTs

Treatment assignment in randomized controlled trials provide a perfect source of instruments. It is fully randomized (complete or conditional), but sometimes people do not fully comply. 

Example: Lottery to head start.

- Head start is an early childhood education (before primary school) program in the U.S.
- In the early 2000s, the U.S. government lauched a lottery program providing access to head-start to eligible families.
- If a family won the lottery, their kid can attend the program.
- It turns out there's no full compliance.
- Can use lottery as an instrument to study the effect of early childhood education on people's later development.
- Kline and Walters (2016), "Evaluating public programs with close substitute: the case of head start", *Quarterly Journal of Economics*.

### 4.2 Natural/Quasi-Experiment

Sometime nature or policy does experiment for you.

- Natural experiment: You cannot control natural phenomenon like earthquakes, typhoons, or even rainfalls. They may affect the treatment but are uncorrelated with all other factors that affect the outcome.
- Quasi-experiment: Policies are usually made without considering people’s idiosyncratic heterogeneity ($U_{i}$), so can be viewed as exogenous.

Example of natural experiment: Rainfall, poverty, and crime.

- Hypothesis: Low-income may incites violence, so usually we can see the regional income level is negatively correlated with the crime rate.
- Regressing crime rate on regional income is problematic. Maybe a high crime rate also causes low-income. 
  - Reverse causality.
- Some economists then look at agriculturally-dependent regions, and use rainfall level as the instrument.
  - Miguel, Staynath and Sergenti (2004), "Economic shocks and civil conflict: an instrumental variables approach", *Journal of Political Economy*.
- Exclusion: Rainfall may not affect crime rate through other channels.
  - Counterexample: Sarsons (2015), "Rainfall and conflict: a cautionary tale", *Journal of Development Economics*, finds that crime rate is still highly correlated with rainfall level in regions whose income is not sensitive to rainfall (e.g., downstream of irrigation dams).

Example of quasi experiment: Maternal care and birthweight.

- Hypothesis: More maternal care (healthcare service for mother-to-be during pregnancy) may causally contribute to high birthweight of newborns. (Birthweight is an important indicator for babies' future cognitive and noncognitive development.)
- Regressing birthweight on maternal care is problematic because women who seek better maternal care may have healthier lifestyles, better economic status, etc. 
  - Self-selection.
- The U.S. government reduced cost/provided free access to maternal care service to the poor. In the late 80s and early 90s, they abruptly lowered the threshold in the definition of being “poor”.
- Exclusion: The policy has nothing to do with any other factors that affect birthweight.
  - Questionable: Different states chose different implementation date, so some families moved to get the benefit. Then we have self-selection again.

### 4.3 More Exogenous Variables than the Treatment

When there’s no way to conduct RCT and no way to find a good natural/quasi-experiment (which is very common), this is almost the last option. The quality of such instruments is usually more questionable than the previous ones, but at least they might tell more the truth than regressing on the treatment directly.

Example: Return to Education by Nobel Laureates.

- Hypothesis: more education causally lead to higher income.
- Regressing income on education is problematic because ability is omitted.
- Instrument 1: Quarter of birth.
  - Angrist and Krueger (1995), "Does compulsory school attendance affect schooling and earnings?", *Quarterly Journal of Economics*.
  - Exclusion: Quarter of birth may be uncorrelated with any cognitive/noncognitive ability of children, and may be also uncorrelated with family background.
  - Relevance: If you have to be 6 years old or above to enter elementary school, those who are born in the last quarter are on average 1 year older than others. They may be more likely to attain less education to earn money, support the family, etc. Could be weak.
- Instrument 2: Distance to a college when children are in high-middle school.
  - David Card (1995), "Using geographic variation in college proximity to estimate the return to schooling", in *Aspects of Labour Market Behacio~lr:Essays in Honour of John Vanderkamp*.
  - Exclusion: Whether you live near a college may be uncorrelated with any cognitive/noncognitive ability of children. It may be correlated with family background, but these data are available and can be included as controls.
  - Relevance: Families who choose to live near a college may have a strong emphasis on education.

Some other common sources: 

- More aggregated level variables: Data are individual-level or household level but use neighborhood level variation as instrument.
  - Example: Household consumption of nutrients, instrumented by consumption of nutrients in reference groups.
  - Dubois, Griffith, and Nevo (2014), "Do prices and attributes explain international differences in food purchases?", *American Economic Review*.
- Hausman-type instrument.
  - Example: Endogenous price instrumented by prices of the same product in other markets (districts).
  - Nevo (2001), "Measuring market power in the ready-to-eat cereal industry", *Econometrica*.
