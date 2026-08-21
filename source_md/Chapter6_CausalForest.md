# ECON5280 Chapter 6 Causal Forest

<font size="5">Junlong Feng</font>

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/junlong-feng/econ5280/main?filepath=Lecture6_Forest.ipynb)

## Outline

* Motivation: ATE only gives us the average but we may want to know heterogeneity.
* CATE: A causal parameter that reflects heterogeneity of ATE in observables.
* Traditional Methods: A brief intro to/review of nonparametrics.
* Random Tree, Causal Tree and Causal Forest: Econometrics learns from and develops ML.
* Applications and Implementation: Everything can be done in R by some simple commands.

## 1. CATE and Heterogeneity

In the last chapter, we learned ATE which is defined as
$$
ATE\equiv \mathbb{E}(Y_{i}(1)-Y_{i}(0))
$$
for binary treatment. With an i.i.d. random sample, ATE is the same for everyone $i$. This is not always desirable since it conceals all heterogeneity; people are different, so exposure to the same treatment may lead to different causal effect.

To capture heterogeneity, the individual treatment effect, $ITE_{i}\equiv Y_{i}(1)-Y_{i}(0)$, of course does the job, but it has two drawbacks: 

- Impossible to estimate due to *the fundamental problem of causal inference*.
- Even if we can estimate it, it does not give us any insights: We would know that the effect of attending to UST on one's future salary could be different between Alice and Bob if $ITE_{Alice}\neq ITE_{Bob}$, but we don't know what leads to the difference: there are so many things not in common for Alice and Bob but the difference in ITE fails to find the link.

In this chapter, we study a new causal parameter which resolves all the previously mentioned issues. Suppose we have a bunch of variables $W_{i}$, we define (still assuming a binary treatment) the **conditional average treatment effect** as follows:
$$
CATE(W_{i})\equiv \mathbb{E}(Y_{i}(1)-Y_{i}(0)|W_{i}).
$$
Now by construction, 

- CATE contains richer information that ATE: for each value $w$ of $W_{i}$, we have a different value of CATE. Moreover, once we know the entire CATE function, we can backout ATE by the law of iterated expectations:
  $$
  ATE=\mathbb{E}(CATE(W_{i})).
  $$

- The richer information contained in CATE is more structured than ITE. Let $W_{i}$ be gender, equal to 1 if $i$ is female. Then by comparing $CATE(1)$ and $CATE(0)$, we know how gender affects the economic return to attending UST, which has more economic implications than simply known that the economic return to attending UST differs between Alice and Bob as ITE does.

### 1.1 Unconfoundedness and Identification of CATE

Similar to ATE, CATE is defined in terms of the potential outcomes, so we cannot directly estimate it. We need to first establish identification to link it to the observables. 

In this chapter we introduce a weaker assumption than **complete randomization** in Chapter 5. It does not require RCT and can be satisfied even in observational data.

**Unconfoundedness**: $D_{i}\perp (Y_{i}(1),Y_{i}(0))|W_{i}$. 

Unconfoundedness says that the treatment only needs to be independent of the potential outcomes conditional on $W_{i}$. To see this is weaker than complete randomization,

- Under complete randomization, $D_{i}\perp (Y_{i}(1),Y_{i}(0),W_{i})$, which implies "unconfoundedness" and $D_{i}\perp W_{i}$. However, here we allow $D_{i}$ to be correlated with $W_{i}$.
- Think about the model $Y_{i}(1)=g_{1}(W_{i},U_{i})$ and $Y_{i}(0)=g_{0}(W_{i},U_{i})$. Complete randomization requires that $D_{i}\perp (U_{i},W_{i})$ but unconfoundedness allows $D_{i}$ to be **correlated with both $W_{i}$ and $U_{i}$$**; it only requires that $D_{i}\perp U_{i}$ conditional on $W_{i}$.

When does unconfoundedness hold? Unconfoundedness can hold in both experimental data and observation data.

- First, since complete randomization implies unconfoundedness mathematically, unconfoundedness holds in RCT data if the treatment is fully randomly assigned.
- In some RCT, randomization is conditional on $W_{i}$ by construction. For instance, suppose I want to study the causal effect of attending tutorials on final grade. Now I offer a tutorial to randomly chosen students from a university. Randomization is essentially made by flipping a coin: If you get a head, you must attend the tutorial. If it's a tail, then you cannot participate. Suppose from an observational study, there are more female students attending tutorials than male students. To make my experiment sample match this proportion so that my result would be more reality-relevant, I will prepare two coins with different head-probability. It has a higher probability of getting a head by flipping coin 1 than coin 2. Then I'll use coin 1 for female students and coin 2 for male students so that more female students will eventually receive the treatment. Let participation of the tutorial be $D_{i}$ where $D_{i}=1$ if $i$ participate. Let $W_{i}=1$ if $i$ is female. Then $D_{i}$ is fully randomized conditional on either $W_{i}=1$ or $W_{i}=0$, so **unconfoundedness** holds. However, $D_{i}$ is correlated with $W_{i}$ by construction, **violating complete randomization**.
- In observational data, unconfoundedness is sometimes called **selection on observables**. It means that $D_{i}$ is not randomized by an experimenter but self-selected by $i$. When $i$ is making the selection, unconfoundedness says her choice is only based on $W_{i}$, having nothing to do with $U_{i}$. For instance, let $W_{i}$ contain family background, undergraduate college, major and GPA, gender, age and MBTI. Let $U_{i}$ contain other earning abilities not captured by $W_{i}$. If we can observe $W_{i}$ in a dataset and if we believe whether coming to UST or not only depends on $W_{i}$, that is, once those factors in $W_{i}$ is fixed, whether coming to UST or not is just a random action, then unconfoundedness holds.
  - A direct implication is that, in observational data, we usually need a relatively big set of variables in $W_{i}$ to justify unconfoundedness because we worry about omitting a variable which the decision $D_{i}$ is based on.

Now given unconfoundedness, we are ready to derive identification of CATE.

**Theorem**. Under unconfoundedness, $CATE(W_{i})=\mathbb{E}(Y_{i}|D_{i}=1,W_{i})-\mathbb{E}(Y_{i}|D_{i}=0,W_{i})$.

*Proof*. By unconfoundedness,
$$
\begin{align*}
CATE(W_{i})\equiv& \mathbb{E}(Y_{i}(1)|W_{i})-\mathbb{E}(Y_{i}(0)|W_{i})\\
\overset{\text{unconfoundedness}}{=}&\mathbb{E}(Y_{i}(1)|D_{i}=1,W_{i})-\mathbb{E}(Y_{i}(0)|D_{i}=0,W_{i})\\
=&\mathbb{E}(Y_{i}|D_{i}=1,W_{i})-\mathbb{E}(Y_{i}|D_{i}=0,W_{i}).
\end{align*}
$$
Now we can think about how to estimate CATE by estimating $\mathbb{E}(Y_{i}|D_{i}=1,W_{i})-\mathbb{E}(Y_{i}|D_{i}=0,W_{i})$. At the first glance, it's just the difference of two conditional expectations. We can approximate them, in principle, by sample averages using the MM estimation approach. However, this turns out much harder than it appears.

- Estimating a conditional expectation when the conditioning variable is discrete is easy. 
  - For instance, let $W_{i}\in\{0,1\}$. Then $\mathbb{E}(Y_{i}|D_{i}=1,W_{i}=1)$ can be estimated by taking the average of $Y_{i}$s for the **subsample** where $D_{i}=1$ and $W_{i}=1$.
- Estimating a conditional expectation when the conditioning variable is continuous is difficult. 
  - When the conditioning variable $W$ is continuous, $W=w$ has probability 0 for any $w$.
  - Impossible (with probability one) to find a **subgroup** $\{i\}$ whose $W_{i}=w$.
  - For instance, suppose $W$ is annual income. Suppose I ask you to estimate $\mathbb{E}(Y|W=10125)$, i.e., the expected $Y$ for someone whose annual income is exactly equal to 10125. It is highly possible that in your data set there does not exist any individual whose income is equal to this number.

## 2. Traditional Methods

Now let's forget about treatment effects for a moment. The conditional expectation is a function of the conditioning variables. For instance, $\mathbb{E}(Y|D=1,W)$ changes when $W$ changes. So if we recover the entire functions: $\mathbb{E}(Y|D=1,W=w)$ and $\mathbb{E}(Y|D=0,W=w)$ for all possible $w$ that $W$ can take, then we are done. For notational simplicity, let
$$
f_{d}(w)\equiv \mathbb{E}(Y|D=d,W=w).
$$
 Our goal is to estimate function $f_{d}(w)$ for different $d$s.

In statistics, generally, there are two approaches to estimate a function. Let's call them **local methods** and **gloabl methods**. 

**Local methods**. Local methods view the function to be estimated as a collection of values: Knowing the whole function $f_{d}(\cdot)$ is equivalent to knowing $f_{d}(w)$ for every value $w$ in the domain. So we can estimate one function value at a time. For an arbitrary fixed $w$, local methods estimate $f_{d}(w)$ using $Y$s around $W=w$. For instance, suppose I assme funtion $f_{d}(w)$ is continuous in $w$, then if I go arbitrarily close to $w$, the function value will also be arbitrarily close to $f_{d}(w)$.

- Suppose $W$ is a scalar and is continuous.
- Suppose we want to know $f_{d}(10)$.
- Ideally we would like to find all $Y_{i}$s whose $D_{i}=d$ and $W_{i}=10$ and then take average.
- But this may not be possible. So instead, we find all $i$s whose $D_{i}=d$ and $W_{i}\in (10-h,10+h)$, and take average of their $Y_{i}$s.
- Of course there is bias because, say, $f_{d}(9.8)$, is not equal to $f_{d}(10)$ even though they are close. 
- Bias is larger if $h$ is larger. Variance is smaller is $h$ is larger because the effective sample size increases.
- If I send $h\to 0$ when $n\to\infty $, the bias should be gone asymptotically and I get a consistent estimate.

Key logic: When $w$ changes a little, $f_{d}(w)$ changes a little. So one needs to find observations whose $W$'s values are close to the target $w$.

Key challenge: **Curse of dimensionality**.

**Curse of dimensionality**. When $W$ is no longer a scalar, containing $p>1$ variables, it's hard to find *nearby* observations:

- Suppose you have $n$ data points, i.e., $n$ individuals.
- Suppose all variables in $W$ take values on $[0,1]$.
- When $p=1$, it's like throwing $n$ points randomly into a unit interval. Points are close to each other.
- When $p=2$, throw $n$ points into a unit square. When $p=3$, a cube. And so on...
- On average, your $n$ data points lie farther away from each other as $p$ gets larger, creating the so-called **sparsity**.
- For instance, suppose $W$ only contains income. Then for $W=10125$, you only need to find people whose income is around 10125. But if $W$ contains income and experience, then for instance, for $W=(10125,10)$, you need to find people whose income is around 10125 **with** around 10 year experience. Fewer people will satisfy both criteria.

Theoretically, now your **effective sample size** is no longer $n$, but (of the order of) $nh^{p}$ where $p$ is the dimensionality of $W_{i}$. Therefore, the standard error of your estimator is large (recall that the se is proportional to $\sqrt{sample\ size}$). Meanwhile, since consistency is driven by the fact that the estimator's variance is shrinking to 0, such an estimator's rate of convergence to the true function value is also slow.

**Global methods**. Global methods, instead of treating a function as a collection of values, treat the function as an entity, a member contained in some functional space $\mathscr{F}$. Think about the finite dimensional Euclidean space $\mathbb{R}^{k}$. There exists $k$ elements in it forming a **basis** such that any $k$ vector can be written as a linear combination of the them. Now suppose $\mathscr{F}$ has similar properties, such as the Hilbert space, then there could also be a elements in it, i.e., functions in it, forming a basis such that my desired function is a linear combination of the basis functions. The problem is, for $\mathbb{R}^{k}$, the number of elements in the basis is $k$. Now $\mathscr{F}$ is infinite-dimensional (why?), so the basis contains infinite functions. Let $f^{(1)},...$ be a basis. Then there exists constants $a_{1},...$ such that
$$
f_{d}=\sum_{j=1}^{\infty}a_{j}f^{(j)}
$$
Unlike $f_{d}(w)$, the basis $f^{(j)}$ functions are known. Popular choices include B-splines, polynomials and wavelets. As long as we can figure out these $a_{j}$s, we know $f_{d}$. 

 Similar to $h$ in local methods, here we also need to do approximation because we cannot handle infinite sum in practice. We need to choose a cutoff $J$, and approximate $f_{d}$ by $\sum_{j=1}^{J}a_{j}f^{(j)}$.  In theory, we let $J\to\infty$ as $n\to\infty$.

- A larger $J$ yields smaller bias but larger variance. A smaller $J$ yields smaller variance but larger bias.
- Curse of dimensionality still exists: The more variables you have, the more complicated the approximation is (think about polynomials).
- Rate of convergence is again slower than $\sqrt{n}$, slowed down by $J\to\infty$.

Note that throughout, we never considered the value of $f_{d}$ at any specific $w$. Global methods deal with the entire function directly.

**Machine learning** provides new methods to estimate functions, but still, fundamentally, they belong to these two types. We will see important algorithms in each type in this semester: We do tree and forest methods in this chapter, which can be thought of as an adaptive local method. We will do neural networks near the end of the semester, which is a global method.

## 3. Random Tree, Causal Tree, and Causal Forest

For local methods, the key problem behind the curse of dimensionality is that we are obsessed with the idea of finding $i$ whose $W$ is close to $w$. However, this is sufficient but not necessary to approximate $f_{d}(w)$.

- It is possible that $f_{d}(w)=f_{d}(w')$ when $w$ and $w'$ are far away. So always only focusing on points near $w$ loses information.
- How near is near depends on the flatness of $f_{d}(w)$ around $w$. As an extreme case, if $f_{d}$ is a constant function, then we can use the full sample to estimate it.

Therefore, perhaps we can find regions of $w$ on which $f_{d}(w)$ is relatively flat. On each region, we just take the average of $Y_{i}$. These regions may be wide or tight; they are determined by data, not by our ad-hoc bandwidth $h$. This is the essence of the tree methods.

### 3.1 Random Tree

Suppose $p=2$, i.e., there are two variables in $W$: $W_{1}$ and $W_{2}$. Further, imagine $(W_{1},W_{2})$ take values on the unit square $[0,1]\times [0,1]$. At the moment, forget about CATE for simplicity. Let's say our goal is simply to estimate a conditional expectation function $f_{d}(w)\equiv \mathbb{E}(Y_{i}|D_{i}=d,W_{i}=w)$.

- Now forget about the target value $w$ and we'll never try to find values close to $w$ any more.
- Instead, we are going to split the unit square into small rectangles.
- The goal is, in each rectangle, the conditional expectation, or, $f_{d}(w)$, is almost a constant.
- After we finish, we go back and find out in which rectangle our target $w$ lies. Then simply average $Y_{i}$s for the $i$s whose $W$ take value in the rectangle.

Two problems:

1. How do we split the unit square efficiently? There are infinite ways to split it.
2. How do we know whether $f_{d}(w)$ changes a lot in a rectangle or not? We need to estimate $f_{d}(w)$ by taking average of the $Y$s in a rectangle, but by construction its estimate does not change at all for all $w$ in the rectangle.

**The random tree algorithm** (a skecth).

1. Split the sample into two halves by $W_{1}<t_{1}$ and $W_{1}>t_{1}$. Calculate the average of $Y_{i}$s in these two subsamples. Denote them by $\bar{Y}_{t1}^{(1)}$ and $\bar{Y}_{t1}^{(2)}$. Calculate $\Delta_{1}(t_{1})\equiv (\bar{Y}_{t1}^{(1)}-\bar{Y}_{t1}^{(2)})^{2}$. 
2. Try all possible $t_{1}\in (0,1)$. Find the one that yields the largest $\Delta_{1}(t_{1})$. Call it $t_{1}^{*}$. Mathematically, $t_{1}^{*}=\arg\max_{t_{1}}\Delta_{1}(t_{1})$.
3. Split the **original** sample into two halves by $W_{2}<t_{2}$ and $W_{2}>t_{2}$. Repeat Steps 1 and 2 for all possible $t_{2}\in (0,1)$ and find $t_{2}^{*}$.
4. Compare $\Delta_{1}(t_{1}^{*})$ and $\Delta_{2}(t_{2}^{*})$. We choose **the regressor**, 1 or 2, and its corresponding $t^{*}$ which yields the largest $\Delta(t^{*})$ as our first-round split rule.
5. Now we have two subsamples. In each of them, repeat Steps 1-4.
6. Repeat the splitting until the pre-set number of rounds splitting and/or the size of the final rectangles are reached. The final rectangles are called nodes or **leaves**.

If you learned random tree from ML, you may remember one needs a testing set. This algorithm is simpler by taking advantage of the splitting criterion we use. The idea is from, e.g., Wager and Athey (2018, JASA) and Athey, Tibshirani and Wager (2019, Annals of Stats). Feel free to explore the difference if you’re interested but this is not required by this course.

This algorithm solves the two questions we raised earlier:

1. How do we split the unit square efficiently?

   - No matter how large $p$ is, each time we only focus on one of them. Each splitting trial is a simple one-dimensional problem. Very easy to handle; the complexity of the problem is linear in $p$.
   - The above algorithm is for $p=2$. When $p>2$, only needs to add similar steps as step 3 for every covariate. Then in step 4, compare $\Delta_{1}(t_{1}^{∗})$, $\Delta_{2}(t_{2}^{∗})$, ..., and $\Delta_{p}(t_{p}^{∗})$, and find the largest one.
   - This **recursive splitting** algorithm is so simple that it can be represented by a binary tree shown below. 

2. How do we know whether $f_{d}(w)$ changes a lot in a leaf or not? 

   - We actually do not look at $f_{d}(W)$ i.e. $E(Y|D=d,W_{1},W_{2})$ on one leaf. We compare the difference between the $f_{d}(W)$s at different leaves. We find the split that yields that largest difference (largest contrast, as some authors prefer). 
   - In this way, the points on one leaf A have relatively similar function value $f_{d}(\cdot)$ compared with the points in another leaf B. Otherwise, those points would have been classfied onto leaf B.
   - When the leaves are small enough, we are confident that function $f_{d}(\cdot)$ is almost constant on it, leading to small enough bias when averaging $Y_{i}$s on it.

<img src="/Users/junlong/Documents/HKUST Teaching/Fall 2022/ECON5280/illuTree.png" alt="illuTree" style="zoom: 50%;" />

Figure from *the Elements of Statistical Learning* by Friedman, Tibshirani and Hastie (Springer, 2009), p.306.

- Top left: A generic partition. Hard to construct. Hard to interpret.
- Top right: A random tree. (I know it doesn’t look like a tree!) Easy to interpret. Hard to draw when $p>2$.
- Bottom left: Standard way to represent a random tree. (Here we go!) Easy to draw for an arbitrarily large $p$.
- Bottom right: Estimated conditional expectation on each leaf. The vertical axis is Y. 

### 3.2 Causal Tree

So far we used tree to estimate conditional expectations. However, our ultimate goals are i) to estimate causal effect and ii) conduct statistical inference. These are two new questions to answer. 

1. How do we estimate CATE?
2. How do we establish statistical properties like consistency and asymptotic distribution?

#### 3.2.1 Response to Q1: Changing the Splitting Critirion

In the tree algorithm, we treated $\mathbb{E}(Y|D=d,W=w)$ as a function of $w$, i.e., $f_{d}(w)$. We set splitting criterion by comparing $\hat{f}_{d}$ for $W_{j}>t_{j}$ and $W_{j}<t_{j}$, and set the split as the $j$ and $t_{j}$ that yields the largest contrast.

Now since our goal is to estimate $\mathbb{E}(Y|D=1,W=w)-\mathbb{E}(Y|D=0,W=w)$, we can simply treat this difference as the function we care about:

- Let $f(w)\equiv \mathbb{E}(Y|D=1,W=w)-\mathbb{E}(Y|D=0,W=w)$.
- Revise Step 1 in the algorithm in [Section 3.1](#3.1 Random Tree) as follows:
  1. Split the sample into two halves by $W_{1}<t_{1}$ and $W_{1}>t_{1}$. In each subsample, calculate the average of $Y$ for those observations with $D=1$ minus the average of $Y$ for those observations with $D=0$. Denote them by $\bar{Y}_{1-0,t1}^{(1)}$ and $\bar{Y}_{1-0,t1}^{(2)}$. Calculate $\Delta_{1}(t_{1})\equiv (\bar{Y}_{1-0,t1}^{(1)}-\bar{Y}_{1-0,t1}^{(2)})^{2}$. 
- Then for all the following steps, revise $\Delta$ in the same way.

In this way, we can guarantee that on each leaf, i.e., the final rectangle, the CATE does not change much in $W$, so averaging the $Y$s on it does not create much bias.

#### 3.2.2 Response to Q2: Honesty

The second problem is much harder to solve. 

- Stefan Wager and Susan Athey in a series of papers solve the problem by introducing a notion called “honesty”.
- Recall in the algorithm, $Y$ is used for two purposes: a) determine the covariate and the split in each round, and b) after the leaves are formed, compute the final estimates of CATE.
- Such dependence causes theoretical challenges for statistical properties.
  - For instance, consistency could be derived by applying WLLN on each leaf: at the end of the day, we just do sample average on each leaf.
  - WLLN needs i.i.d.
  - However, these leafs are formed by comparing $Y$s. So each leaf is a function of $Y$. Conditional on the leaf, $Y$ are no longer i.i.d.

To resolve this issue, Wager and Athey propose an extra step before Step 1 in the algorithm in [Section 3.1](#3.1 Random Tree), called *honest splitting*:

**The causal tree algorithm** (a skecth).

- Step 0 (Honest splitting). Randomly split the data set into two halves by $i$ (**NOT BY $W$**!)
  - Each subsample has around $n/2$ observations. Call them the **training data** and **estimation data** respectively.
- Step 1-6. Grow a random tree **using the training data ONLY** following the algorithm in [Section 3.1](#3.1 Random Tree) with the modified $\Delta$ in [Section 3.2.1](#3.2.1 Changing the Splitting Critirion). 
- Given the leaves, for $w$ of interest, estimate $CATE(w)$ by first check which leaf $w$ falls into, and then calculate $\bar{Y}_{D=1}-\bar{Y}_{D=0}$ **using the estimation data**.

*Honesty* means that we use independent subsamples to grow a tree and to estimate, respectively. Then the bias caused by “Y is correlated to Y” can be avoided.

- Wager and Athey show that under some regularity conditions, the estimated CATE is unbiased, consistent and asymptotically normal (**but, like most nonparametric estimators, at a slower rate of convergence than $1/\sqrt{n}$**). Inference can be easily done by standard t-test, p-value and confidence interval.

### 3.3 Causal Forest

A causal tree has two drawbacks:

1. It’s sensitive to noises.
2. Some data (the estimation data) are never used to grow a tree, while some (the training data) are never used for estimation. Information loss.

**Causal Forest**: Grow many trees and take the average.

**The causal forest algorithm** (a skecth).

1. Randomly draw $s$ observations from the full sample of $n$ observations.
2. Grow a causal tree using this subsample with the algorthm in [Section 3.2.2](#3.2.2 Honesty).
   - The steps are exactly what we described earlier. Recall the key steps are: splitting this subsample into two halves (honesty), recursively splitting the training data to grow a tree, and use the estimation data to estimate the CATE.
3. Repeat Steps 1 and 2 for $B$ times. Choose $B$ as large as possible.
4. You end up with B causal trees. They are called a **causal forest**. For any given $w$, each tree $b$ yields an estimate of $CATE^{(b)}(w)$. Take average of theses $B$ estimates: $\sum_{b=1}^{B}CATE^{(b)}(w)/B$.

The causal forest resolves the two drawbacks of causal trees by:

1. It’s sensitive to noises.
   - Averaging makes the noise in each tree only have $1/B$ share of impact.
2. Some data (the estimation data) are never used to grow a tree, while some (the training data) are never used for estimation. Information loss.
   - Training data in one tree may be used as estimation data in another tree. When B is large enough, all data are likely to be used to grow trees and to estimate.

Wager and Athey show that under certain regularity conditions: CATE estimated by causal forest is consistent and asymptotically normal (**again, with a slower rate of convergence than $\sqrt{n}$ **).

#### 3.3.1 Implementation Details

Pay attention to the following three aspects in implementation.

**Covariates subsetting**. 

Recall that when growing a tree, in each round in the recursive splitting, one needs to try every possible split for every covariate. 

- This can be slow when $p$ is large.
- In practice, for each round in the recursive splitting, one only need to try out a random subset of covariates. 
- This subset is randomly chosen.
- A rule of thumb is randomly choosing $\min⁡\{\sqrt{p}+20,p\}$ covariates.

**Minimum leaf size**.

Leaf size is like $h$ in traditional nonparametrics. A large leaf usually leads to a large bias (because extrapolated too much). A small leaf may result in large variance (because data points are too few).

- We can control the minimum size of the leaf: when it is reached, recursive splitting ends.
- This size can be cross-validated.

**Imbalance of a split**.

- When splitting a parent node, the size of each child node is not allowed to be too different.
- Meanwhile, the number of treated and untreated observations in a child node can not be too different either.
- The level of imbalance can be controlled in the algorithm.

### 3.4 Comparison with OLS

We saw in Chapter 5 that OLS regression of $Y$ on $D$ or on $(D,W)$ leads to a consistent estimate of the ATE regardless of the nonlinearity of the true model. However, this is no longer true when it comes to CATE. The reason is that CATE is a function of $W$ instead of a value, so misspecificaltion about how $W$ enters the model may lead to a different function of $W$ by construction. Moreover, if $D$ is not independent of $W$, allowed by unconfoundedness, then leaving some $W$ in the error term may lead to violation of the zero correlation assumption between $D$ and the error.

Causal forest overcomes these drawbacks because i) it does not impose any restrictions on the functional form and ii) everything works under conditional randomization of $D$ whereas $W$ can be arbitrarily correlated with the unobservable(s) or potential outcomes.

## 4. Implementation and Application 

You can imagine there is a wide applications of causal forest in economics. Whenever you have a D$ satisfying the unconfoundedness assumption and you want to know the heterogeneous effects, you can replace linear models (or, any parametric nonlinear models) with it and obtain everything you want. (Everything includes ATE; we'll study one of the most powerful estimators for ATE under unconfoundedness in the next chapter.)

Implementation is straightfoward as an R package is ready at. Check out the tutorial webpage https://grf-labs.github.io/grf/index.html. 

Here's an example.

- I generate a data set with $p=10$ and $n=2000$.
- The $2000\times 10$ matrix of $W$ are randomly drawn from normal.
- Treatment $D$ is conditionally randomized: It's binary and equal to 1 with probability 0.4 if $X_{1}<0$ and with probability 0.6 if $X_{1}>0$.
- The model for $Y$ is nonlinear: $Y=\max\{W_{1},0\}\times D+W_{2}+\min\{W_{3},0\}+\varepsilon$.
  - Can you calculate the true CATE by hand? Try.

```R
library(grf)
library(DiagrammeR)
### Generate data. You do not need to run these if a dataset is already given. 
n <- 2000      
p <- 10
W <- matrix(rnorm(n * p), n, p)
D <- rbinom(n, 1, 0.4 + 0.2 * (W[, 1] > 0))
Y <- pmax(W[, 1], 0) * D + W[, 2] + pmin(W[, 3], 0) + rnorm(n)
data=data.frame(Y,D,W)

## Generate values of interest for w
W.test <- matrix(0, 101, p)
W.test[, 1] <- seq(-2, 2, length.out = 101)

### Build a causal forest.
tau.forest <- causal_forest(W, Y, D) # Put regressors
tree <- get_tree(tau.forest, 1)
plot(tree) # visualize the first tree in the forest

### Estimate CATE at specified w values using the forest
tau.hat <- predict(tau.forest, W.test) 
plot(W.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions, 0, 2), 
     xlab = "W1", ylab = "CATE", type = "l") # estimated CATE
lines(W.test[, 1], pmax(0, W.test[, 1]), col = 2, lty = 2) # true CATE

# Estimate treatment effects with confidence interval
tau.forest <- causal_forest(W, Y, D, num.trees = 4000)
tau.hat <- predict(tau.forest, W.test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)
plot(W.test[, 1], tau.hat$predictions, ylim =
       range(tau.hat$predictions + 1.96 * sigma.hat, 
             tau.hat$predictions - 1.96 * sigma.hat, 0, 2), 
     xlab = "W1", ylab = "CATE", type = "l")
lines(W.test[, 1], tau.hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(W.test[, 1], tau.hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)
lines(W.test[, 1], pmax(0, W.test[, 1]), col = 2, lty = 1)

```

Next, we estimate CATE using the data set in Chapter 5. 

```R
data <- read.csv("welfare-small.csv")
n <- nrow(data)
Y=as.matrix(data$Y,n,1)
D=as.matrix(data$D,n,1)
W=cbind(data$age,data$educ)


## Generate values of interest for w
W.test <- cbind(30,seq(min(data$educ),max(data$educ),by=1))

### Build a causal forest.
tau.forest <- causal_forest(W, Y, D) # Put regressors
tree <- get_tree(tau.forest, 1)
plot(tree) # visualize the first tree in the forest

# Estimate treatment effects with confidence interval
tau.forest <- causal_forest(W, Y, D, num.trees = 4000)
tau.hat <- predict(tau.forest, W.test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)
plot(W.test[, 2], tau.hat$predictions, ylim =
       range(tau.hat$predictions + 1.96 * sigma.hat, 
             tau.hat$predictions - 1.96 * sigma.hat), 
     xlab = "educ", ylab = "CATE", type = "l")
lines(W.test[, 2], tau.hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(W.test[, 2], tau.hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)

```













 

