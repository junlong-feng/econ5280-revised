# ECON5280 Chapter 11 Difference-in-Differences and Event Studies

<font size="5">Junlong Feng</font>

## Outline

* Motivation: Policies are often adopted at different times in different places. Can we use these timing differences to learn causal effects?
* The Two-Group and Two-Period Difference-in-Differences Design: Identification by parallel trends.
* Staggered Treatment Adoption: Define group-time average treatment effects when different groups are treated at different times.
* What Can Go Wrong with Two-Way Fixed Effects: Already-treated groups may become invalid controls.
* Modern Difference-in-Differences and Event Studies: Estimate transparent causal effects first and aggregate them second.
* Credibility and Inference: What event-study graphs can and cannot tell us.
* Implementation in R.

## 1. Motivation

Many important treatments in economics are policies rather than randomized experiments. A policy may be introduced in one place first and in another place later. For example:

* Some cities introduce a housing subsidy while other cities do not.
* Some provinces raise the minimum wage earlier than other provinces.
* Some firms adopt a new technology while other firms have not adopted it yet.

Suppose we want to study the effect of a housing subsidy on rent. A simple comparison of rent between cities with and without the subsidy is unlikely to be causal.

* Cities adopting the subsidy may already have higher rent before adoption.
* These permanent differences do not disappear simply because we have panel data.

Another simple idea is to compare rent in an adopting city before and after adoption. This comparison is also unlikely to be causal.

* Rent may change over time even without the subsidy because of inflation, interest rates, or a common housing boom.

Difference-in-differences (DID) combines the two comparisons:

1. Calculate the change over time for the treated group.
2. Calculate the change over time for an untreated comparison group.
3. Subtract the second change from the first.

The permanent difference between the two groups is removed by the first difference. The common change over time is removed by the second difference. The remaining difference can be interpreted as a treatment effect under a **parallel trends assumption**.

Before introducing notation, consider a numerical example. Suppose average monthly rent is:

| Group | Before | After | Change |
|---|---:|---:|---:|
| Cities adopting the subsidy | 120 | 134 | 14 |
| Comparison cities | 100 | 106 | 6 |

The adopting cities have higher rent in both periods, so comparing levels is misleading. Rent rises by 14 in the adopting cities but by 6 in the comparison cities. The DID estimate is therefore
$$
(134-120)-(106-100)=14-6=8.
$$
Under parallel trends, 6 is our estimate of how much rent would have risen in the adopting cities without the subsidy. The remaining 8 is the estimated effect on the adopting cities.

## 2. Two Groups and Two Periods

We begin with the simplest DID design. There are two time periods, $t=0$ and $t=1$, and two groups.

* $G_i=1$: Unit $i$ belongs to the treated group. This group is untreated at $t=0$ and treated at $t=1$.
* $G_i=0$: Unit $i$ belongs to the control group. This group is untreated in both periods.
* $D_{it}=G_i\times 1(t=1)$ is the realized treatment status.
* $Y_{it}(1)$ is the potential outcome under the regime "treated beginning in period 1," and $Y_{it}(0)$ is the potential outcome under the regime "never treated."
* The observed outcome is
  $$
  Y_{it}=G_iY_{it}(1)+(1-G_i)Y_{it}(0).
  $$

The main parameter is the average treatment effect on the treated in the post-treatment period:
$$
ATT\equiv \mathbb{E}\left[Y_{i1}(1)-Y_{i1}(0)\mid G_i=1\right].
$$

Why do we focus on ATT instead of ATE? Policy adoption is not randomized. DID most directly recovers the average effect for the group that actually adopts the policy. The effect for a group that never adopts is generally a different parameter and may not be identified by the same design.

### 2.1 No Anticipation

Before the policy starts, units should not already be affected by it:
$$
Y_{i0}(1)=Y_{i0}(0).
$$

In the simple two-period design, this means the observed pre-treatment outcome is untreated for both groups:
$$
Y_{i0}=Y_{i0}(0).
$$

No anticipation can fail. For instance, firms may start hiring before a previously announced subsidy formally begins. In that case, the last nominally pre-treatment period may already be affected by treatment.

### 2.2 Parallel Trends

The key identifying assumption is
$$
\mathbb{E}\left[Y_{i1}(0)-Y_{i0}(0)\mid G_i=1\right]
=
\mathbb{E}\left[Y_{i1}(0)-Y_{i0}(0)\mid G_i=0\right].
$$

This is called **parallel trends**. It says that, without treatment, the average outcome of the treated group would have changed by the same amount as the average outcome of the control group.

Some important remarks:

* Parallel trends concerns **changes**, not levels. The treated and control groups may have very different average outcomes before treatment.
* Parallel trends does **not** say that treatment is independent of the potential outcomes. Selection into treatment based on permanent differences is allowed.
* Parallel trends concerns the untreated potential outcome $Y(0)$. After treatment begins, the untreated path of the treated group is unobserved, so the assumption cannot be verified directly.

### 2.3 Identification

Under no anticipation and parallel trends,
$$
\begin{align*}
ATT
=&\ \mathbb{E}\left[Y_{i1}(1)-Y_{i1}(0)\mid G_i=1\right]\\
=&\ \mathbb{E}\left[Y_{i1}-Y_{i0}\mid G_i=1\right]
-\mathbb{E}\left[Y_{i1}(0)-Y_{i0}(0)\mid G_i=1\right]\\
=&\ \mathbb{E}\left[Y_{i1}-Y_{i0}\mid G_i=1\right]
-\mathbb{E}\left[Y_{i1}-Y_{i0}\mid G_i=0\right].
\end{align*}
$$

The last line is the DID identification equation. Every term in it is an expectation of observed outcomes.

The corresponding estimator is
$$
\widehat{ATT}_{DID}
=\left(\bar{Y}_{11}-\bar{Y}_{10}\right)
-\left(\bar{Y}_{01}-\bar{Y}_{00}\right),
$$
where $\bar{Y}_{gt}$ is the sample average outcome for group $g$ in period $t$.

The order of the two differences does not matter:
$$
\underbrace{(\bar{Y}_{11}-\bar{Y}_{10})-(\bar{Y}_{01}-\bar{Y}_{00})}_{\text{change over time, then difference across groups}}
=
\underbrace{(\bar{Y}_{11}-\bar{Y}_{01})-(\bar{Y}_{10}-\bar{Y}_{00})}_{\text{difference across groups, then change over time}}.
$$

### 2.4 Regression Representation

With two groups and two periods, the DID estimator is also the OLS coefficient $\hat{\tau}$ from
$$
Y_{it}=\beta_0+\beta_1G_i+\beta_2Post_t+\tau(G_i\times Post_t)+\varepsilon_{it},
$$
where $Post_t=1(t=1)$. The interaction $G_i\times Post_t$ is exactly the treatment indicator $D_{it}$.

For panel data, the same coefficient can be obtained from
$$
Y_{it}=\alpha_i+\lambda_t+\tau D_{it}+\varepsilon_{it},
$$
where $\alpha_i$ is a unit fixed effect and $\lambda_t$ is a time fixed effect. This is called a **two-way fixed effects (TWFE)** regression.

In this simple two-group and two-period setting, TWFE is just another way to calculate DID. The difficulty discussed later does not arise yet.

### 2.5 Panel Data and Repeated Cross Sections

There are two common data structures.

* **Panel data**: The same units are observed before and after treatment. We can form $Y_{i1}-Y_{i0}$ for every unit.
* **Repeated cross sections**: Different random samples are observed in each period. We cannot form an individual change, but we can still compare the four group-time averages.

Both designs can identify ATT under appropriate sampling and parallel-trends assumptions. Panel data are often more efficient because permanent unit-level noise disappears after differencing.

## 3. Staggered Treatment Adoption

The simple design is useful, but many empirical applications have more than two periods and different groups adopt treatment at different times.

Let $t=1,\ldots,T$ and define

* $G_i=g$: Unit $i$ is first treated in period $g$.
* $G_i=\infty$: Unit $i$ is never treated. In computer code, never-treated units are often recorded as $G_i=0$.
* $D_{it}=1(t\geq G_i)$ for eventually treated units and $D_{it}=0$ for never-treated units.

We focus on **staggered adoption**: once treatment begins, it remains on. Hence,
$$
D_{it}=1\implies D_{i,t+1}=1.
$$

This setup does not cover a treatment that repeatedly switches on and off without additional assumptions.

### 3.1 Potential Outcomes and Group-Time ATT

Let $Y_{it}(g)$ be the potential outcome in period $t$ if unit $i$ is first treated in period $g$, and let $Y_{it}(0)$ be the outcome if unit $i$ is never treated.

No anticipation requires
$$
Y_{it}(g)=Y_{it}(0),\qquad t<g.
$$

Treatment effects may differ across adoption groups and calendar periods. The natural building block is
$$
ATT(g,t)\equiv
\mathbb{E}\left[Y_{it}(g)-Y_{it}(0)\mid G_i=g\right],\qquad t\geq g.
$$

For example:

* $ATT(4,4)$ is the average effect in the adoption period for units first treated in period 4.
* $ATT(4,6)$ is the average effect two periods after adoption for the same group.
* $ATT(6,6)$ concerns a different group in the same calendar period as the previous example.

There is no reason these three effects must be equal. A policy may take time to work, and early adopters may differ from late adopters.

### 3.2 Choosing a Valid Control Group

To estimate $ATT(g,t)$, we need a nonempty group that is untreated in both the baseline period $g-1$ and period $t$. Two common choices are:

* **Never-treated controls**: Compare group $g$ with units that never receive treatment.
* **Not-yet-treated controls**: Compare group $g$ with units whose treatment begins after $t$, together with any never-treated units.

The second choice uses more observations but requires us to believe that later adopters provide a valid counterfactual trend for earlier adopters. The choice of control group is therefore an identifying decision, not merely a software option.

For never-treated controls, the relevant parallel-trends assumption is
$$
\mathbb{E}\left[Y_{it}(0)-Y_{i,g-1}(0)\mid G_i=g\right]
=
\mathbb{E}\left[Y_{it}(0)-Y_{i,g-1}(0)\mid G_i=\infty\right],
$$
with $\Pr(G_i=\infty)>0$. Under no anticipation, this implies

$$
ATT(g,t)=
\mathbb{E}\left[Y_{it}-Y_{i,g-1}\mid G_i=g\right]
-
\mathbb{E}\left[Y_{it}-Y_{i,g-1}\mid G_i=\infty\right].
$$

For not-yet-treated controls, replace the right side of the parallel-trends assumption by
$$
\mathbb{E}\left[Y_{it}(0)-Y_{i,g-1}(0)\mid G_i>t\right]
$$
and require $\Pr(G_i>t)>0$. This gives
$$
ATT(g,t)=
\mathbb{E}\left[Y_{it}-Y_{i,g-1}\mid G_i=g\right]
-
\mathbb{E}\left[Y_{it}-Y_{i,g-1}\mid G_i>t\right],
$$
where $G_i>t$ includes never-treated units by convention.

Each $ATT(g,t)$ is therefore constructed from a clean two-group and two-period DID comparison. Already-treated observations are not used as untreated controls.

### 3.3 Conditional Parallel Trends

Unconditional parallel trends may be implausible when adoption groups have different covariate distributions. Let $W_i$ contain pre-treatment covariates. A conditional version is
$$
\begin{align*}
&\mathbb{E}\left[Y_{it}(0)-Y_{i,g-1}(0)\mid G_i=g,W_i\right]\\
&\qquad=
\mathbb{E}\left[Y_{it}(0)-Y_{i,g-1}(0)\mid G_i>t,W_i\right].
\end{align*}
$$

We also need **overlap**: for every covariate value represented in adoption group $g$, there must be not-yet-treated controls with that covariate value. Formally,
$$
\Pr(G_i>t\mid W_i=w)>0
\quad\text{for every }w\in\operatorname{supp}(W_i\mid G_i=g).
$$
When using never-treated controls, replace $G_i>t$ by $G_i=\infty$. Without overlap, the control group's conditional trend is unavailable for some treated units.

Under conditional parallel trends and overlap,
$$
\begin{align*}
ATT(g,t)=\mathbb{E}\Big[&
\mathbb{E}(Y_{it}-Y_{i,g-1}\mid G_i=g,W_i)\\
&-\mathbb{E}(Y_{it}-Y_{i,g-1}\mid G_i>t,W_i)
\ \Big|\ G_i=g\Big].
\end{align*}
$$

The inner expectations compare units with the same $W$. The outer expectation averages the conditional effects over the covariate distribution of the treated group $G_i=g$.

This connects DID to earlier chapters:

* Outcome regression estimates the conditional change for the control group.
* Inverse propensity-score weighting reweights the controls to resemble adoption group $g$.
* A doubly robust DID estimator combines the two approaches and remains consistent if one of the two nuisance models is correctly specified, under the relevant regularity conditions.

Covariates should generally be measured before treatment and should not themselves be caused by treatment. Controlling for a post-treatment variable can destroy the causal interpretation.

## 4. Why Conventional TWFE Can Fail

With many periods, a familiar regression is
$$
Y_{it}=\alpha_i+\lambda_t+\beta D_{it}+u_{it}.
$$

It is tempting to interpret $\beta$ as an average treatment effect. This interpretation is valid in the canonical two-group and two-period design and under some strong homogeneous-effect models. It is not generally valid with staggered adoption and heterogeneous effects.

### 4.1 The Bad-Control Problem

Consider three groups:

* Group E is treated early.
* Group L is treated late.
* Group N is never treated.

When group E is newly treated, groups L and N are untreated and can both be reasonable controls. When group L is newly treated, however, group E is already treated. A conventional TWFE regression may use all three kinds of comparisons:

1. Newly treated versus never treated: potentially valid.
2. Newly treated versus not yet treated: potentially valid.
3. Newly treated versus already treated: generally invalid as an untreated comparison.

The third comparison subtracts a change that may contain the earlier group's treatment-effect dynamics. It does not reveal the change that the newly treated group would have experienced without treatment.

Consequently, the TWFE coefficient can combine causal effects using difficult-to-interpret or negative implicit weights. It is useful to distinguish two statements. The Goodman--Bacon decomposition assigns nonnegative weights to its component two-by-two DID comparisons. But when TWFE is rewritten as a weighted average of the underlying cohort-time treatment effects, some implicit weights can be negative. Moreover, a comparison using an already-treated group can subtract changes in that group's treatment effect. In extreme cases, every underlying treatment effect can be positive while the TWFE coefficient is negative.

### 4.2 Conventional Event-Study Regressions

Applied researchers often estimate
$$
Y_{it}=\alpha_i+\lambda_t+
\sum_{e\neq -1}\mu_e 1(t-G_i=e)+u_{it},
$$
where $e=t-G_i$ is **event time**:

* $e=-2$: two periods before treatment.
* $e=-1$: the omitted reference period.
* $e=0$: the treatment-adoption period.
* $e=2$: two periods after treatment.

With simultaneous treatment timing, this is a natural extension of DID. Under staggered adoption and heterogeneous effects, however, a coefficient $\mu_e$ may be contaminated by treatment effects at other event times. Even a lead coefficient can be nonzero because of post-treatment effect heterogeneity rather than a genuine pre-trend.

Therefore, simply adding leads and lags to TWFE does not solve the bad-control problem.

### 4.3 When Is TWFE Still Useful?

The lesson is not that fixed effects are always wrong.

* In the two-group and two-period design, TWFE equals the DID estimator.
* With one treated cohort sharing a common treatment date and a valid never-treated control group, conventional event-study regressions do not use already-treated groups as controls. If every unit is treated on the same date, event-time indicators are absorbed by calendar-time fixed effects and the effects are not identified.
* Under sufficiently strong treatment-effect homogeneity assumptions, the static TWFE coefficient can recover the common effect.
* Unit and time fixed effects remain useful ingredients inside modern estimators.

The lesson is that, under staggered adoption, the researcher must define the causal parameter and the valid controls before running a regression.

## 5. Modern DID and Event-Study Estimators

Modern methods follow a simple principle:

> Estimate transparent cohort-time causal effects using untreated controls, and only then aggregate them into the desired summary.

### 5.1 Group-Time Effects

The Callaway--Sant'Anna approach directly estimates the collection of $ATT(g,t)$ values. Depending on the assumptions and available covariates, each clean comparison can be estimated by:

* difference-in-means of outcome changes,
* outcome regression,
* inverse probability weighting, or
* a doubly robust DID estimator.

The collection of $ATT(g,t)$ values is informative but can be large. Aggregation should match the empirical question.

### 5.2 Overall, Group, and Calendar-Time Effects

An overall effect averages post-treatment $ATT(g,t)$ values using explicit nonnegative weights. One possible summary is
$$
\theta^{overall}=\sum_{g}\sum_{t\geq g}\omega_{g,t}ATT(g,t),
\qquad \omega_{g,t}\geq 0,\quad \sum_{g,t}\omega_{g,t}=1.
$$

Other aggregations answer different questions:

* **Group aggregation**: How large is the average effect for each adoption cohort?
* **Calendar-time aggregation**: What is the average effect in each calendar period?
* **Overall aggregation**: What is a single summary effect among treated observations?

These are different parameters. Software cannot decide which one is economically relevant.

In the `did` package used below, `type = "simple"` weights each post-treatment $ATT(g,t)$ in proportion to cohort size. Early cohorts then receive more total weight because they contribute more post-treatment cells. `type = "group"` first averages post-treatment effects within each cohort and then averages across cohorts. `type = "dynamic"` averages the contributing cohorts at each event time using cohort-size weights. These summaries need not agree because they answer different questions.

### 5.3 Dynamic Effects and Event Studies

For event time $e=t-g$, define
$$
\theta(e)=
\sum_{g:g+e\leq T}
\Pr(G_i=g\mid G_i+e\leq T)ATT(g,g+e).
$$

For $e\geq 0$, $\theta(e)$ is the average effect after exactly $e$ periods of exposure among cohorts observed at that event time. Plotting $\theta(e)$ against $e$ gives a heterogeneity-robust event study.

Be careful when comparing distant event times. Early adopters can be observed many periods after treatment, whereas late adopters cannot. The set of cohorts contributing to $\theta(e)$ can therefore change with $e$. A changing event-study curve may reflect both:

* genuine treatment-effect dynamics, and
* a changing composition of adoption cohorts.

One solution is to report a balanced event window, using only cohorts observed for every event time in the chosen window. This improves comparability but discards data.

### 5.4 Other Modern Approaches

Several modern estimators implement the same broad principle in different ways.

* **Sun--Abraham interaction-weighted estimator**: Estimate cohort-by-event-time interactions and aggregate them, avoiding contamination from other event times.
* **Borusyak--Jaravel--Spiess imputation estimator**: Estimate the untreated outcome model using untreated observations, impute untreated counterfactuals for treated observations, and average the resulting treatment-effect estimates.
* **Group-time DID estimators**: Construct clean two-by-two comparisons for each $(g,t)$ and then aggregate.

The estimators can differ in their target parameters, comparison groups, efficiency, covariate adjustment, and maintained models. They should not be treated as interchangeable commands. The research design should determine the estimator.

## 6. Credibility, Pre-Trends, and Inference

### 6.1 What a Pre-Trend Plot Can Tell Us

Before treatment, a modern event study can report **placebo** or **pseudo-ATT** estimates. Large systematic differences in pre-treatment changes are evidence against the proposed parallel-trends design.

However, estimates close to zero do not prove parallel trends.

* Pre-treatment estimates may be noisy and have low power.
* Untreated paths may be parallel before treatment but diverge afterward.
* Selecting a specification only because its pre-trend test is insignificant changes the statistical behavior of the final estimates.
* The identifying assumption concerns an unobserved post-treatment counterfactual, not merely the observed pre-treatment path.

Pre-trends are therefore a diagnostic, not a certificate of validity.

### 6.2 Simultaneous Confidence Bands

An event study reports many coefficients. Pointwise 95% confidence intervals cover each coefficient separately with probability approximately 95%, but they do not provide 95% coverage for the entire curve.

For a joint statement such as “all pre-treatment effects are zero” or for visual inspection of the complete path, simultaneous confidence bands are preferable. Modern DID packages can construct uniform bands using multiplier bootstrap methods.

### 6.3 Clustering

Outcomes for the same unit are correlated over time. If treatment is assigned at a higher level, outcomes can also be correlated within that assignment cluster.

* With state-level policies, standard errors are usually clustered at the state level rather than the individual level.
* With few treated or control clusters, conventional cluster-robust approximations may be unreliable.

The clustering level should follow the source of treatment assignment and dependence, not the number of rows in the data set.

### 6.4 A Design Checklist

Before estimating a staggered DID, answer the following questions.

1. What is the unit and what is the treatment-adoption date?
2. Is treatment absorbing, or can it turn off?
3. Could units anticipate treatment?
4. Which observations are valid controls for each treated cohort and period?
5. Is parallel trends more plausible unconditionally or conditional on pre-treatment covariates?
6. What causal parameter is of interest: group-time, dynamic, calendar-time, or overall ATT?
7. Does the composition of cohorts change across the event-time window?
8. At what level should inference be clustered?
9. Are other policies or shocks introduced at the same time?
10. Is the conclusion sensitive to the outcome scale, sample window, control group, and anticipation window?

## 7. Implementation in R

We illustrate staggered DID using simulated panel data. Some units are first treated in periods 4, 6, or 8; others are never treated. Treatment effects are dynamic and differ across adoption cohorts. This is exactly the situation in which conventional TWFE can be difficult to interpret.

### 7.1 Generate the Data

```R
set.seed(5280)

N <- 1200                         # number of units
T <- 10                           # number of periods

# G=0 denotes never treated. Other values are first-treatment periods.
units <- data.frame(
  id = 1:N,
  G = sample(c(0, 4, 6, 8), N, replace = TRUE,
             prob = c(0.30, 0.25, 0.25, 0.20)),
  alpha = rnorm(N),
  W = rnorm(N)
)

# Long panel: one row for every unit-period pair.
data <- merge(expand.grid(id = 1:N, t = 1:T), units,
              by = "id", sort = TRUE)

data$D <- as.integer(data$G > 0 & data$t >= data$G)
data$event_time <- ifelse(data$G > 0, data$t - data$G, NA)

# Early and late adopters have different treatment-effect paths.
data$cohort_scale <- ifelse(data$G == 4, 1.50,
                     ifelse(data$G == 6, 1.00,
                     ifelse(data$G == 8, 0.60, 0)))

data$tau <- ifelse(data$D == 1,
                   data$cohort_scale * (1 + 0.5 * data$event_time), 0)

# Untreated outcomes have unit effects and common time effects.
data$lambda <- 0.25 * data$t + 0.03 * data$t^2
data$Y0 <- data$alpha + data$lambda + 0.5 * data$W + rnorm(nrow(data))
data$Y <- data$Y0 + data$tau
```

In this DGP, untreated outcomes satisfy parallel trends. The treatment effect is positive but differs across cohorts and grows with exposure time.

### 7.2 Calculate One Group-Time DID by Hand

The following function estimates $ATT(g,t)$ using period $g-1$ as the baseline and units not yet treated by $t$ as controls.

```R
simple_att_gt <- function(data, g, t) {
  stopifnot(t >= g)
  base <- g - 1

  now <- data[data$t == t, c("id", "G", "Y")]
  names(now)[names(now) == "Y"] <- "Y_now"

  before <- data[data$t == base, c("id", "Y")]
  names(before)[names(before) == "Y"] <- "Y_before"

  change <- merge(now, before, by = "id")
  change$dY <- change$Y_now - change$Y_before

  treated <- change$G == g
  controls <- change$G == 0 | change$G > t

  mean(change$dY[treated]) - mean(change$dY[controls])
}

# Effect for cohort 4 in period 6: two periods after adoption.
simple_att_gt(data, g = 4, t = 6)

# True effect in this simulation, shown only because the DGP is known.
mean(data$tau[data$G == 4 & data$t == 6])
```

The two numbers should be close in a large sample. With real data, the second number is an unobserved counterfactual parameter and cannot be calculated.

### 7.3 Estimate All Group-Time Effects

The `did` package implements Callaway and Sant'Anna's group-time approach. Its treatment-timing variable must equal the first treated period for treated units and 0 for never-treated units.

```R
library(did)

att <- att_gt(
  yname = "Y",
  tname = "t",
  idname = "id",
  gname = "G",
  xformla = ~ 1,
  data = data,
  panel = TRUE,
  control_group = "notyettreated",
  est_method = "dr",
  base_period = "universal",
  clustervars = "id",
  bstrap = TRUE,
  cband = TRUE
)

summary(att)
ggdid(att)
```

Because this example has no covariate adjustment (`xformla = ~ 1`), outcome-regression, IPW, and doubly robust implementations reduce to the same unconditional two-by-two DID point estimates. With pre-treatment covariates, one can replace `~ 1` by, for example, `~ W`.

The option `base_period = "universal"` uses period $g-1$ as the common baseline and normalizes event time $-1$ to zero, matching the event-study notation above. The package default, `base_period = "varying"`, leaves post-treatment estimates unchanged but uses adjacent-period comparisons for pre-treatment pseudo-ATTs.

The options `bstrap = TRUE` and `cband = TRUE` request multiplier-bootstrap simultaneous confidence bands. In an application where treatment is assigned at a level above `id`, add the appropriate assignment-level cluster where supported (for example, `clustervars = c("id", "state")`). The package permits at most two clustering variables and requires one of them to be `idname` for panel data.

### 7.4 Aggregate the Effects

```R
# A single weighted average of post-treatment group-time effects.
overall <- aggte(att, type = "simple")
summary(overall)

# Average effects by adoption group.
by_group <- aggte(att, type = "group")
summary(by_group)

# Heterogeneity-robust event study.
dynamic <- aggte(att, type = "dynamic", min_e = -3, max_e = 4)
summary(dynamic)
ggdid(dynamic)

# Hold cohort composition fixed through post-treatment event time 3.
dynamic_balanced <- aggte(att, type = "dynamic",
                          min_e = -3, max_e = 3,
                          balance_e = 3)
summary(dynamic_balanced)
ggdid(dynamic_balanced)
```

The `simple`, `group`, and `dynamic` outputs answer different questions. They need not be numerically equal even when every estimate is correctly calculated.

### 7.5 Compare with Conventional and Interaction-Weighted Event Studies

The following optional code compares a conventional TWFE event study with the Sun--Abraham interaction-weighted estimator. The `fixest` function `sunab()` creates cohort-by-event-time interactions and aggregates them appropriately.

```R
library(fixest)

data$ever_treated <- as.integer(data$G > 0)
data$event_twfe <- ifelse(data$G > 0, data$t - data$G, 0)

# Conventional TWFE event study.
twfe_event <- feols(
  Y ~ i(event_twfe, ever_treated, ref = -1) | id + t,
  data = data,
  cluster = ~ id
)

# fixest treats a cohort outside the observed time range as never treated.
data$G_sunab <- ifelse(data$G == 0, 1000, data$G)

# Sun--Abraham interaction-weighted event study.
sa_event <- feols(
  Y ~ sunab(G_sunab, t, ref.p = -1) | id + t,
  data = data,
  cluster = ~ id
)

iplot(list("Conventional TWFE" = twfe_event,
           "Sun-Abraham" = sa_event),
      ref.line = -1,
      main = "Event-study estimates")
```

In this DGP, treatment effects differ across cohorts and exposure lengths. Compare the two curves and relate any difference to the controls and implicit aggregation used by each estimator.

The `did` plots above use multiplier-bootstrap simultaneous confidence bands because `bstrap = TRUE` and `cband = TRUE`. The intervals displayed by this `fixest::iplot()` comparison are ordinary pointwise intervals, so do not interpret them as simultaneous bands for the whole curve.

## Summary

* DID identifies causal effects by comparing outcome changes, not outcome levels.
* The key assumptions are no anticipation and parallel trends for untreated potential outcomes.
* With two groups and two periods, the DID estimator equals the interaction coefficient in a regression and the treatment coefficient in TWFE.
* Under staggered adoption, $ATT(g,t)$ is the natural causal building block.
* Valid controls are never-treated or, under an appropriate assumption, not-yet-treated units. Already-treated units are generally not valid untreated controls.
* Conventional static and dynamic TWFE regressions can be difficult to interpret when treatment effects are heterogeneous.
* Modern methods estimate clean group-time or cohort-event-time effects first and aggregate them second.
* An event-study graph is useful for displaying dynamics and diagnosing pre-trends, but insignificant pre-treatment estimates do not prove parallel trends.
* Aggregation, event-window balance, anticipation, covariates, and clustering are substantive design choices.

## Further Readings

* Callaway, B. and Sant'Anna, P. H. C. (2021), [Difference-in-Differences with Multiple Time Periods](https://doi.org/10.1016/j.jeconom.2020.12.001), *Journal of Econometrics*.
* Sun, L. and Abraham, S. (2021), [Estimating Dynamic Treatment Effects in Event Studies with Heterogeneous Treatment Effects](https://doi.org/10.1016/j.jeconom.2020.09.006), *Journal of Econometrics*.
* Goodman-Bacon, A. (2021), [Difference-in-Differences with Variation in Treatment Timing](https://doi.org/10.1016/j.jeconom.2021.03.014), *Journal of Econometrics*.
* Borusyak, K., Jaravel, X., and Spiess, J. (2024), [Revisiting Event-Study Designs: Robust and Efficient Estimation](https://doi.org/10.1093/restud/rdae007), *Review of Economic Studies*.
* Roth, J., Sant'Anna, P. H. C., Bilinski, A., and Poe, J. (2023), [What's Trending in Difference-in-Differences? A Synthesis of the Recent Econometrics Literature](https://doi.org/10.1016/j.jeconom.2023.03.008), *Journal of Econometrics*.
* Sant'Anna, P. H. C. and Zhao, J. (2020), [Doubly Robust Difference-in-Differences Estimators](https://doi.org/10.1016/j.jeconom.2020.06.003), *Journal of Econometrics*.
* Rambachan, A. and Roth, J. (2023), [A More Credible Approach to Parallel Trends](https://doi.org/10.1093/restud/rdad018), *Review of Economic Studies*.
* R implementation: [`did` package documentation](https://bcallaway11.github.io/did/) and [`fixest::sunab()` documentation](https://lrberge.github.io/fixest/reference/sunab.html).
