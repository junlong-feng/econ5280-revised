# Chapter 5 Causal Inference and Randomized Controlled Trials (RCT)

<font size="5">Junlong Feng</font>

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/junlong-feng/econ5280/main?filepath=Lecture5_RCT.ipynb)

## Outline

* Motivation: Econ cares about causality, but what is causality? How is it possible to make causal inference?
* Potential Outcomes: The key elements in Rubin's approach to causal inference.
* Average Treatment Effect (ATE): The first causal parameter we study in this semester.
* RCT and Identification of ATE: Understand the most powerful technique to conduct causal inference and then express ATE in terms of expectations of observables.
* Estimator 1: Difference-in-means.
* Estimator 2: Regression adjustment.

## 1. Potential Outcomes

There are more than one way to formally define causality. This course follows the one that is mostly adopted in frequentist statistics and econometrics: Rubin's potential outcome framework (developed in the 1970s).

From now on, let's imagine we have someone called $i$:

* Treatment $D_{i}$: A random variable that indicates the treatment that $i$ gets.
  * e.g. Whether $i$ goes to UST or not./ Whether $i$ studies for a master's degree or not.
  * $D_{i}$ can be either discrete or continuous in general, but we focus on the binary case in this chapter. 
  * $D_{i}$ in upper case is random. The lower case $d$ represents the nonrandom treatment status. For instance, let $d\in\{0,1\}$, and $D_{i}=0$ means $i$ does not go to UST whereas $D_{i}=1$ means $i$ goes to UST.
* Potential outcome $Y_{i}(d)$: A random variable representing the outcome that $i$ would get if she took treatment $d$.
  * e.g. Corresponding to the first example above, $Y_{i}(1)$ may represent $i$'s would-be future income if she went to UST and $Y_{i}(0)$ is then $i$'s would-be income if she did not.
  * Once we write the potential outcome by $Y_{i}(d)$, we have already made assumptions: i) We implicitly assume that other people's treatment status does not affect $i$'s potential outcome and ii) we implicitly assume that everyone can have treament status $d$. These two assumptions are called **stable unit treatment value assumption (SUTVA)** and we'll maintain it implicitly throughout the semester.
* Observed or realized outcome $Y_{i}$: $Y_{i}=Y_{i}(d)$ if and only if $D_{i}=d$. Or equivalently, $Y_{i}=Y_{i}(D_{i})$.
  * The realized outcome is **the potential outcome** associated with **the** treatment status which $i$ indeed receives.
  * e.g. Suppose you are $i$. You have already come to UST to study. So your $d=1$, and your future income will be equal to your potential future income at $d=1$. However, there may exist a *parallel universe* where you didn't come to UST (i.e. your $d=0$) . Then the actual future income of *that you* in the parallel world will be $Y_{i}(0)$.
* Individual causal, or, treatment effect (**ITE**), of receiving treament $d$ compared to $d'$: $Y_{i}(d)-Y_{i}(d')$.

## 2. The Average Treatment Effect (ATE)

Although mathematically well-defined, ITE is almost impossible to estimate. This is due to the so called *fundamental problem of causal inference*:

**Fundamental problem of causal inference**: You only live once so you only know at most one of $Y_{i}(d)$ and $Y_{i}(d')$. Suppose $D_{i}=1$ so you observe $Y_{i}(1)$. To know the effect of $D_{i}$ on $i$, one needs $Y_{i}(0)$, which is completely unknown. How about using someone else' ($j$ for example) outcome to replace the potential outcome of $i$ that I cannot observe? 

* Suppose your friend $j$ didn't come to UST $(d=0)$, and you want to use her income, i.e., $Y_{j}(0)$ as your $Y_{i}(0)$. But this is actually assuming that $Y_{i}(0)=Y_{j}(0)$, which is highly unlikely: You are assuming $i$ and $j$ would have **exactly** the same income if they both did not attend UST. We know that income is potentially affected by a lot of factors, say, IQ, ability, age, family background, and perhaps, unfortunately, gender. Even if all these inputs are miraculously equal for $i$ and $j$, the income production functions of the two individuals also need to be equal. 

So, let us give up ITE and pursue a more approachable causal parameter, the average treatment effect **ATE**: 
$$
ATE(d,d')\equiv \mathbb{E}(Y_{i}(d))-\mathbb{E}(Y_{i}(d'))
$$

* Note ATE has NO subscript because we are going to work with an i.i.d. sample, so $i$'s expectation is equal to everyone else and thus subscript $i$ can be dropped.
* As mentioned, we focus on the case where $D$ is binary in this chapter. Let $d=1$ and $d'=0$, the only two values that $D$ can take on. So from now on, we denote $ATE(d,d')$ in short by $ATE$.

## 3. RCT and Identification of ATE as the Difference in Means

RCTs are often considered as the gold standard for causal inference. The key condition that RCT buys us is the **complete randomization**:
$$
D_{i}\perp(Y_{i}(1),Y_{i}(0)).
$$
Let's see what an RCT is and why complete randomization holds via an example:

- Suppose a new drug is developed. Scientists want to understand whether the drug is effective or not. 
  - Recruit a representative group of people (experiment subjects).
  - Randomly split the group into two subgroups: Treatment group takes the drug, denoted by $D=1$, whereas the control group gets a placebo, denoted by $D=0$.
  - Let $Y_{i}(d)$ be whether $i$ fully recovers 3 days after taking the drug $(d=1)$/placebo $(d=0)$. Then ATE: $\mathbb{E}(Y_{i}(1))-\mathbb{E}(Y_{i}(0))$.
  - By construction, whether individual $i$ takes the drug or not is not selected by herself but instead randomly assigned to her, which has nothing to do with her preference/health condition/family backgroupd, etc. Mathematically, random variables $Y_{i}(1)$ and $Y_{i}(0)$ are thus fully independent of $D_{i}$, i.e., $D_{i} \perp (Y_{i}(1),Y_{i}(0))$. This is the most important property RCT gives us. We will show that it plays a central role in identifying and estimating ATE.
  - **Be careful**: If it's not an RCT, for instance, there was no experiment done and we survey people who took the drug when they were ill and people who did not. By asking whether they recovered 3 days after the onset of the illness, $Y$, we can have a similar data set. However, in this case, **it is more likely that** $D_{i}\not\perp (Y_{i}(1),Y_{i}(0))$ because whether $i$ took the drug or not may be correlated with her wealth, her lifestyle, her baseline health condition, etc. All these factors may impact $Y_{i}(d)$ as well, resulting in nonzero correlation between $D_{i}$ and $Y_{i}(d)$s.

We just saw that RCT buys us complete randomization. Now let's see what complete randomization buys us. It turns out that under complete randomization, we can express ATE, a parameter defined in terms of the **potential outcomes**, in terms of expectations that are only related to **observed variables**:
$$
\begin{align*}
ATE\equiv &\mathbb{E}(Y_{i}(1))-\mathbb{E}(Y_{i}(0))\\
=&\mathbb{E}(Y_{i}(1)|D_{i}=1)-\mathbb{E}(Y_{i}(0)|D_{i}=0)\ \ \ \ (\text{by }D_{i} \perp (Y_{i}(1),Y_{i}(0)))\\
=&\mathbb{E}(Y_{i}|D_{i}=1)-\mathbb{E}(Y_{i}|D_{i}=0)\ \ \ \ (\text{by }Y=Y(d)\text{ if and only if }D=d).
\end{align*}
$$
The last equality achieves **identification of ATE**: ATE is defined in terms of the potential outcomes, but in the end we express it in terms of the moments of the observed variables: $Y$ and $D$.

Now that ATE is identified, it is possible to estimate it using a data set. We introduce two different estimators. They are both unbiased, consistent and asymptotically normal. The tradeoff is that, the first estimator is slightly easier to compute whereas the second has a smaller standard error, that is, more efficient. 

## 4. Estimator 1: Difference-in-means

Given that $ATE=\mathbb{E}(Y_{i}|D_{i}=1)-\mathbb{E}(Y_{i}|D_{i}=0)$, it is very easy to come up with the following estimator. For simplicity, denote $ATE$ by $\tau$. Suppose there are $n_{1}$ individuals who receive $D=1$ and $n_{0}$ receiving $D=0$. Then define
$$
\hat{\tau}_{DM}\equiv \frac{1}{n_{1}}\sum_{D_{i}=1}Y_{i}-\frac{1}{n_{0}}\sum_{D_{i}=0}Y_{i},
$$
where $\sum_{D_{i}=d}$ is a short-hand notation for $\sum_{i\in\{i:D_{i}=d\}},d=0,1$. Each of the two subsample averages are estimating $\mathbb{E}(Y_{i}|D_{i}=1)$ and $\mathbb{E}(Y_{i}|D_{i}=0)$, respectively. We can of course directly derive its statistical properties such as unbiasedness. However, since we have learned regression, I now show you that $\hat{\tau}_{DM}$ actually coincides with the OLS estimator of a linear regression.

Consider $\mathbb{E}(Y|D)$. It is a random variable whose randomness only comes from $D$. Since $D$ is binary, we have
$$
\begin{align*}
\mathbb{E}(Y|D)=&D\mathbb{E}(Y|D=1)+(1-D)\mathbb{E}(Y|D=0)\\
=&\mathbb{E}(Y|D=0)+D\left[\mathbb{E}(Y|D=1)-\mathbb{E}(Y|D=0)\right]\\
=&\mathbb{E}(Y|D=0)+\tau D,
\end{align*}
$$
where the last equality comes from identification, built on our complete randomization condition $D\perp (Y(1),Y(0))$. Now we **define** $\varepsilon\equiv Y-\mathbb{E}(Y|D)$, so by construction we have
$$
\mathbb{E}(\varepsilon|D)=0.
$$
Denote the nonrandom constant $\mathbb{E}(Y|D=0)$ by $\delta_{0}$, so we automatically have
$$
Y_{i}=\delta_{0}+\tau D_{i}+\varepsilon_{i},\mathbb{E}(\varepsilon_{i}|D_{i})=0.
$$
**Remember we DID NOT say that this is the TRUE MODEL of $Y$**. $Y$ can follow any complicated model. But no matter what the true model for $Y$ is, we just proved that **there exists** a random variable $\varepsilon_{i}$ satisfying $\mathbb{E}(\varepsilon_{i}|D_{i})=0$ such that the ATE $\tau$ is the coefficient on $D_{i}$ in linear regression (5). In this proof, the key step relies on **complete randomization**.

We now show that our $\hat{\tau}_{DM}$ coincides with the OLS estimator of $\tau$ in regression (5). Let $X=(1,D)'$. Then $X_{i}X'_{i}=\begin{pmatrix}1&D_{i}\\D_{i}&D^{2}_{i}\end{pmatrix}$. Denote $\beta=(\delta_{0},\tau)'$. By the OLS formula:
$$
\begin{align*}
\hat{\beta}=&\left(\sum_{i}X_{i}X_{i}'\right)^{-1}\left(\sum_{i}X_{i}Y_{i}\right)\\
=&\frac{1}{n\sum_{i}D_{i}^{2}-\left(\sum_{i}D_{i}\right)^{2}}\begin{pmatrix}\sum_{i}D_{i}^{2}\sum_{i}Y_{i}-\sum_{i}D_{i}\sum_{i}D_{i}Y_{i}\\
n\sum_{i}D_{i}Y_{i}-\sum_{i}D_{i}\sum_{i}Y_{i}\end{pmatrix}\\
=&\frac{1}{n_{1}n_{0}}\begin{pmatrix}n_{1}\sum_{i}(1-D_{i})Y_{i}\\
n\sum_{i}D_{i}Y_{i}-n_{1}\sum_{i}Y_{i}\end{pmatrix}\\
=&\begin{pmatrix}\frac{1}{n_{0}}\sum_{D_{i}=0}Y_{i}\\
\frac{1}{n_{1}}\sum_{D_{i}=1}Y_{i}-\frac{1}{n_{0}}\sum_{D_{i}=0}Y_{i}\end{pmatrix},
\end{align*}
$$
so $\hat{\beta}_{2}=\hat{\tau}_{MD}$.

Hence, we do not need to prove that $\hat{\tau}_{DM}$ has properties such as unbiasedness; all the nice properties of $\hat{\tau}_{MD}$ immediately follows from the general theory of OLS estimators because although we constructed it in a different way, we just proved that $\hat{\tau}_{MD}$ is exactly equal to an OLS estimator derived from regression (5). In particular, we have
$$
\begin{align*}
\mathbb{E}(\hat{\tau}_{DM})=\tau,\ \ \ \ &(\text{Unbiasedness})\\
\hat{\tau}_{DM}\overset{p}{\to}\tau,\ \ \ \ &(\text{Consistency})\\
\sqrt{n}(\hat{\tau}_{DM}-\tau)\overset{d}{\to}N\left(0,V_{DM}\right),\ \ \ \ &(\text{Asymptotic Normality})
\end{align*}
$$
where $V_{DM}\equiv\frac{Var(Y|D=1)}{\Pr(D=1)}+\frac{Var(Y|D=0)}{\Pr(D=0)}$. (Can you prove why the asymptotic variance has this expression?)

Finally, for inference, we need to figure out the expression for standard error (se). We can estimate $V_{DM}$ simply by the sample analogue. R can report it directly. Once we have the se, we can do t-test, p-value, and confidence intervals immediately.

**Summary**. There are two ways to calculate the difference-in-means estimator. One can either do it by taking the sample average difference as in (3), or by running a regression of $Y$ on $(1,D)$ (i.e., (5)) and taking the coefficient on $D$. The two approaches yield exactly the same estimates, and there statistical properties are thus identical. 

### 4.1 An Empirical Example

This example is taken from course materials from Susan Athey and Stefan Wager [*Machine Learning-based Causal Inference Tutorial*](https://bookdown.org/stanfordgsbsilab/ml-ci-tutorial/). The data set is a simplified version of General Social Surveys (GSS) ([Smith, 2016](https://gss.norc.org/Documents/reports/project-reports/GSSProject%20report32.pdf)). The setting is an RCT. Individuals were asked about their thoughts on government spending on the social safety net. Social safety nets are non-contributory transfer programs, targeted to the poor or those vulnerable to shocks. Examples include cash transfers, food distribution, ect. The treatment is the wording of the question: about half of the individuals were asked if they thought government spends too much on "welfare" ($D_{i}=1$), while the remaining half was asked if they thought government spends too much on  "assistance to the poor" ($D_{i}=0$). The outcome is also binary, with $Y_{i}=1$ corresponding to "yes".

```R
library(lmtest)
library(sandwich)
# Read in data
data <- read.csv("welfare-small.csv")
n <- nrow(data)

# The difference-in-means estimator
Y=data$Y
D=data$D
tau_dm=mean(Y[D==1])-mean(Y[D==0])

# The least square estimator
model=lm(Y~D)
coeftest(model,vcov=vcovHC(model,type='HC2'))
```

This example shows that (i) the OLS estimator yields exactly the same estimate as the DM estimator, and (ii) wording of the question leads to a significant difference in people's interpretation of government spending.

## 5. Estimator 2: Regression Adjustment

We alreay have an  estimator for ATE under an RCT setting, but, can we do better?

The DM estimator is already unbiased and consistent, so if we can do any better, the only possibility is to propose an estimator that (i) is still unbiased and consistent, and (ii) has a smaller asymptotic variance. With a smaller variance, the standard error will be smaller, leading to a larger t-value, a smaller p-value, and a shorter confidence interval with everything else equal, making statistical inference sharper.

In any branch of data sciences, there are two ways in general to reduce variance of an estimator: (i) introduce more structures and (ii) utilize more information from data. This is because, a large variance means that the data are noisy so the estimator based on them is inaccurate. With more structures about the model or more variation from the data, we are injecting more information into the estimation process, so accuracy of the estimator is improved, resulting in a smaller variance for the estimator.

- Example. If I ask you to make a guess (i.e., an estimator) about tomorrow's weather in Tokyo, without any information, your answer is more or less a random guess. You may give 10 totally different answers if you are asked 10 times, which means your guess has a large variance. Now if I give you a mathematical model which can calculate Tokyo's weather based on Hong Kong's (**more structure**) or I give you the daily data of the weather condition in Tokyo in the past 30 days (**more information from data**), you may give me a more stable guess which is not likely to change no matter how many times you're asked. This means the variance of your estimator is smaller.
- Note that here we are only talking about variance. It does not have any implications on the bias. For instance, the mathematical model about the weather in Tokyo and Hong Kong can be completely wrong, and thus it will lead you to a fully biased estimator with zero variance. 

Now coming back to our problem, having more structures means that we are imposing more assumptions on the functional forms about, for instance, how $Y(d)$ is generated. This is not desirable because models may be wrong. In this course, we try our best to stay as general as possible by imposing only the minimum assumptions on the model. So let's see how **extra information from data** can help us reducing the variance of our estimator.

Let's think about our data set. A data set may contain many variables, but we are only using $Y$ and $D$ at the moment. Suppose there are other variables $W$ which may affect $Y(d)$, then ignoring them is wasting information. The question is, how can one utilize $W$ to enhance efficiency?

Recall that by defining $\varepsilon_{i}=Y_{i}-\mathbb{E}(Y_{i}|D_{i})$, we have $Y_{i}=\delta_{0}+\tau D_{i}+\varepsilon_{i}$. Now let's **force** $W$ to enter the model linearly. Denote $\eta_{i}=\varepsilon_{i}-W_{i}'\gamma$ for some vector $\gamma$. Then we have
$$
Y_{i}=\delta_{0}+\tau D_{i}+W'_{i}\gamma+\eta.
$$
Under complete randomization, we have $D_{i}\perp W_{i}$ as well. It's important to note that here we are ***NOT assuming*** that the true model for $Y_{i}$ is linear in $D_{i}$ and $W_{i}$. Everything is just by construction and some simple algebra. Without loss of generality, assume that $\mathbb{E}(W_{i})=0$; if not, $\delta_{0}$ will just absorb $\mathbb{E}(W_{i}')\gamma$, making no difference in $\tau$. Now, we are ready to prove the following claim **no matter what the true model for $Y_{i}$ is**:

**Claim**. Under complete randomization, there **always** exists a **unique** vector $\gamma$ such that $\mathbb{E}(W_{i}\eta_{i})=0$, provided that $\mathbb{E}(W_{i}W'_{i})$ is full-rank.

To prove it, simply substitute the definition of $\eta$,
$$
0=\mathbb{E}(W_{i}\eta_{i})=\mathbb{E}(W_{i}\varepsilon_{i}-W_{i}W'_{i}\gamma)\implies \gamma=\left[\mathbb{E}(W_{i}W'_{i})\right]^{-1}\mathbb{E}(W_{i}\varepsilon_{i}).
$$
Meanwhile, by construction, $\mathbb{E}(\eta_{i})=0$ and $\mathbb{E}(\eta_{i} D_{i})=0$. Therefore, we can consistently the ATE $\tau$ simply by OLS regressing $Y$ onto $(1,D_{i},W_{i})$ no matter whether this is the true model or not! Let's call this estimator of ATE $\hat{\tau}_{OLS}$. We now show that $\hat{\tau}_{OLS}$ has a smaller asymptotic variance than $\hat{\tau}_{DM}$. 

- Let $X^{new}_{i}=(1,D_{i},W'_{i})'$.  By the general theory of OLS, we have
  $$
  \left(\hat{\delta}_{0},\hat{\tau}_{OLS},\hat{\gamma}'\right)'\overset{d}{\to}N\left(0,\left[\mathbb{E}(X^{new}_{i}X^{new'}_{i})\right]^{-1}\mathbb{E}(\eta^{2}_{i}X^{new}_{i}X^{new'}_{i})\left[\mathbb{E}(X^{new}_{i}X^{new'}_{i})\right]^{-1}\right).
  $$

- Let $X_{i}=(1,D_{i})'$. By $\mathbb{E}(W_{i})=0$ and $D_{i}\perp W_{i}$, we have
  $$
  \begin{align*}
  \left[\mathbb{E}(X^{new}_{i}X^{new'}_{i})\right]^{-1}&=\begin{pmatrix}\mathbb{E}(X_{i}X'_{i})&0\\0&\mathbb{E}(W_{i}W'_{i})\end{pmatrix}^{-1}\\
  &=\begin{pmatrix}\left[\mathbb{E}(X_{i}X'_{i})\right]^{-1}&0\\0&\left[\mathbb{E}(W_{i}W'_{i})\right]^{-1}\end{pmatrix}.
  \end{align*}
  $$

- Substituting it into the variance formula in (9), by block-diagonality, we have
  $$
  \left(\hat{\delta}_{0},\hat{\tau}_{OLS}\right)'\overset{d}{\to}N\left(0,\left[\mathbb{E}(X_{i}X'_{i})\right]^{-1}\mathbb{E}(\eta^{2}_{i}X_{i}X'_{i})\left[\mathbb{E}(X_{i}X'_{i})\right]^{-1}\right).
  $$
  This is exactly the asymptotic distribution of the OLS estimator of the following infeasible regression:
  $$
  Y_{i}-W'_{i}\gamma=\delta_{0}+\tau D_{i}+\eta_{i}.
  $$
  It is infeasible because $\gamma$ is unknown so we cannot construct a dependent variable to run this regression. However, it has theoretical value because we only need to study the variance of the OLS estimator of $\tau$ in the infeasbible regression, which is simpler because there 's no $W$ on the right-hand side.

- Let $\tilde{Y}_{i}\equiv Y_{i}-W'_{i}\gamma$. By comparing the infeasible regression with the regression in Section 4, all we have to do is to replace $Y_{i}$ in the formula of $V_{DM}$ by $\tilde{Y}_{i}$ to obtain the asymptotic variance of $\hat{\tau}_{OLS}$, denoted by $V_{OLS}$:
  $$
  V_{OLS}=\frac{Var(\tilde{Y}_{i}|D_{i}=1)}{\Pr(D_{i}=1)}+\frac{Var(\tilde{Y}_{i}|D_{i}=0)}{\Pr(D_{i}=0)}.
  $$

- By $Y_{i}=\delta_{0}+\tau D_{i}+\varepsilon_{i}$, we have $Var(Y|D=d)=Var(\varepsilon|D=d)$. For $Var(\tilde{Y}_{i}|D_{i}=d)$, we have
  $$
  \begin{align*}
  &Var(\tilde{Y}_{i}|D_{i}=d)\\
  =&Var(\varepsilon_{i}-W'_{i}\gamma|D_{i}=d)\\
  =&Var(\varepsilon_{i}|D_{i}=d)+\gamma'Var(W_{i})\gamma-2cov(\varepsilon_{i},W'_{i}|D_{i}=d)\gamma\\
  =&Var(\varepsilon_{i}|D_{i}=d)+\mathbb{E}(\varepsilon_{i} W_{i})'\left[\mathbb{E}(W_{i}W'_{i})\right]^{-1}Var(W_{i})\left[\mathbb{E}(W_{i}W'_{i})\right]^{-1}\mathbb{E}(\varepsilon_{i} W_{i})\\
  &-2cov(\varepsilon_{i},W'_{i}|D_{i}=d)\left[\mathbb{E}(W_{i}W'_{i})\right]^{-1}\mathbb{E}(\varepsilon_{i} W_{i}).\end{align*}
  $$

- By $\mathbb{E}(W_{i})=0$, we have $\mathbb{E}(W_{i}W'_{i})=Var(W_{i})$ and $cov(\varepsilon_{i},W'_{i}|D_{i}=d)=\mathbb{E}(\varepsilon_{i} W'_{i}|D_{i}=d)=\mathbb{E}(\varepsilon_{i} W'_{i})$. 
  
- Therefore,  
  $$
  Var(\tilde{Y}_{i}|D_{i}=d)=Var(Y_{i}|D_{i}=d)-\mathbb{E}(\varepsilon_{i} W_{i})'\left[\mathbb{E}(W_{i}W'_{i})\right]^{-1}\mathbb{E}(\varepsilon_{i} W_{i}).
  $$
  
- Since $[\mathbb{E}(W_{i}W_{i}')]^{-1}$ is p.d., $V_{OLS}\leq V_{DM}$, where equality holds if and only if $\mathbb{E}(\varepsilon_{i}W_{i})=0$, i.e., $\mathbb{E}(Y_{i}W_{i})=0$, that is, $W_{i}$ is not correlated with $Y_{i}$. This is intuitive because in that case, including $W_{i}$ in the regression only introduces useless information.
  
  - In the example about Tokyo's weather, this is like to say that I give you gas price in Sydney for the past 30 days but still ask you to guess the weather in Tokyo tomorrow. Now you do have extra information, but the information is possibly useless for your task. Your guess of the weather of Tokyo tomorrow is still a random guess with large variance.


Before we end this section, let me emphasize again that we **never assumed** a linear model. Our estimator is **AlWAYS CORRECT** no matter whether the relationship between $Y$ and $(D,W)$ is linear or not.

### 5.1 The Empirical Example Revisited 

We still use the GSS data in Section 4.1 to reestimate ATE using our second estimator.

```R
# Regress Y on D and W
model1=lm(Y~D+age+polviews+income+educ+marital+sex,data)
coeftest(model1,vcov=vcovHC(model1,type='HC2'))
```

### 5.2 Monte Carlo Simulations

We use simulated data to illustrate that $\hat{\tau}_{OLS}$ and $\hat{\tau}_{DM}$ are both close to the true parameter but the former has a smaller variance.

```R
library(lmtest)
library(sandwich)
# Fix sample size = 1000
n=1000

# 100 simulation replications
B=100

# for replication
set.seed(5280)

# Store the estimates
tau_dm=0
tau_ols=0

# Loop over all simulated samples
for (b in 1:B){
		w=rnorm(n,2,1) # create W; here the mean is not 0
		e=rnorm(n,0,1) # create the error term
		d=(rnorm(n,0,1)>0) # create a binary treatment
		y=1+d+w^2*d+exp(w)+e # true model. x does not enter the model linearly
		model1=lm(y~d) # dm estimator
		tau_dm[b]=model1$coefficients[2]
		model2=lm(y~d+w) # regression adjustment
		tau_ols[b]=model2$coefficients[2]
}
# examine the mean and variance
mean(tau_dm)
mean(tau_ols)
var(tau_dm)
var(tau_ols)
```

## Summary

In this chapter, we learned the following:

- What is the potential outcome framework and how are causal effects defined based on it?
- RCT and complete randomization.
- Two estimators for ATE under RCT. Both are consistent and asymptotically normal. The second one utilizing more information from data can achieve smaller variance as long as the extra information is relevant. Both estimators can be computed based on regressions. Importantly, we never make any assumptions on the true functional form of $Y(d)$; it can be arbitrarily nonlinear. We also did not assume $W$ to be independent/uncorrelated with $Y(d)$ or whatsoever. The only assumptions we made in the whole chapter are SUTVA, i.i.d. sampling, complete randomization of the treatment $D$, and $\mathbb{E}(WW')$ being full-rank.
  - RCT is a channel to provide complete randomization, but not the only one. Natural/quasi-experiments may also do the job.

- Practical implications of the theories in this chapters are: If you believe your binary $D$ is independent of the potential outcomes, you can simply regress $Y$ on $D$ and a bunch of covariates $W$ linearly without worrying model misspecification. The coefficient on $D$ is **ALWAYS** ATE. 

