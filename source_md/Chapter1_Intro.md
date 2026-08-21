# Chapter 1 Intro to Metrics              

<font size="5">Junlong Feng</font>

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/junlong-feng/econ5280/main?filepath=Lecture1_Intro.ipynb)

## 1. What is Econometrics?

<font size="2">  *This section heavily builds on Chapter 1 of Econometrics by Bruce Hansen (2022). Please refer to the book chapter for more details.*</font>

The term *econometrics* was created by Ragnar Frisch (1895-1973), a Norwegian economist, one of the principal founders of the Econometric Society, first editor of *Econometrica*, and co-winner of the first Nobel Prize in Economics Science in 1969. He defined the term *econometrics* as follows:

> A word of explanation regarding the term econometrics may be in order. Its definition is implied in the statement of the scope of the [Econometric] Society, in Section I of the Constitution, which reads: “The Econometric Society is an international society for the advancement of economic theory in its relation to statistics and mathematics.... Its main object shall be to promote studies that aim at a unification of the theoretical-quantitative and the empirical-quantitative approach to economic problems....”
>
> But there are several aspects of the quantitative approach to economics, and no single one of these aspects, taken by itself, should be confounded with econometrics. Thus, econometrics is by no means the same as economic statistics. Nor is it identical with what we call general economic theory, although a considerable portion of this theory has a defininitely quantitative character. Nor should econometrics be taken as synonomous with the appli- cation of mathematics to economics. Experience has shown that each of these three viewpoints, that of statistics, economic theory, and mathematics, is a necessary, but not by itself a sufficient, condition for a real understanding of the quantitative relations in modern economic life. It is the *unification* of all three that is powerful. And it is this unification that constitutes econometrics.
>
> ​                                                    												-- Ragnar Frisch, *Econometrica*, 1933, 1, pp. 1-2.

From this definition, econometrics is, roughly, economic theory + math + statistics. The definition has evolved over time. Today, one would say econometrics is, roughly, economic models + statistics + economic data.

So, econometrics deals with data. But what's the difference between econometrics and other data-related disciplines, say, statistics (stats) or machine learning (ML)? We'll try to answer this question partly in this lecture and throughout the semester.

## 2. An Econometric Perspective to Describe Data

Suppose we want to understand whether having a master's degree would increase your future salary. How should we proceed?

Intuitively, we might say, well, let's look at other people's stories:

* Alice has a master's degree and found a job with a monthly salary = \$12K.
* Bob found a job right after graduating from college with a monthly salary = \$8K.
* So I believe having a master's degree would increase **my future salary** by \$4K.
* Or, I believe having a master's degree would at least increase **future salary in general**. 

This thought experiment is not reliable. The most salient problems are:

* Alice and Bob might not be comparable: Bob might graduate 20 years before Alice; gender difference; major; ranking of their colleges, etc.
* The sample size is too small: Alice and Bob could be outliers.

Therefore, methods to resolve the above issues include i) for each individual, collect more  information (more **variables**), and ii) collect more individuals' information (increase the **sample size**). When we do have more information from more individuals and more variables, it is important to have a clear and simple way to represent the information. 

### 2.1 Variables

We usually use uppercase latin letters to represent variables. For instance,

* Let $Y$ be monthly salary. 
* To indicate whose monthly salary we are talking about, we can put a subscript on the variable: $Y_{Alice}$ indicates Alice's salary. More generally, $Y_{i}$ indicates some individual $i$'s salary.
* We also want to represent the variables that could potentially **explain** $Y_{i}$: whether individual $i$ has a master degree ($X_{1i}$) and her/his age ($X_{2i}$).  For notation simplicity, we put all of them into a **vector**: $X_{i}\equiv (X_{1i},X_{2i})'$. The dash represents transpose.
  + Conventionally, all vectors in econometrics, if not stated otherwise, are treated as **column vectors**.

### 2.2 Data Set as A Matrix

Now suppose we run a survey with $n=1000$ participants.  We ask them to report their values of the 3 variables: $Y$ and $X$. Then we will have a **sample** with three **variables** and 1000 **observations**. We represent the sample by $\{Y_{i},X_{i}\}_{i=1,\ldots,n}$ or sometimes $\{(Y_{i},X_{i}):i=1,\ldots,n\}$.

To visualize a sample, or, a data set, let's look at the following R example:

```R
set.seed(5280) # Set seed for replication
n=4   # Sample size. Chose a small one to view data easily.
Y=exp(rnorm(n,0,1)) # generate some fake income data.
X1=rbinom(n,1,0.5) # 1=having an MA degree. 0=no MA degree
X2=sample(22:30, n, replace=T) # generate integers between 22 and 30 for age

# Option 1: organize the variables as a matrix.
data1=cbind(Y,X1,X2) 
# Option 2: organize the variables as a data set
data2=data.frame(income=Y,MA=X1,age=X2)  

View(data1)
View(data2)
```

By running the code, we can see that a data set, mathmatically, is equivalent to a matrix. Each observation $i=1,\ldots,n$ is a row and there are in total $n$ rows. Each variable is a column. 

Throughout the semester, we will use vectors and matrices for representation and  for all mathematical derivations. This could be challenging and takes time to get used to. But it's absolutely worth it. Matrices are handy tools without which many simple results require much more effort to get. The next lecture note reviews matrix algebra.

### 2.3 A Probabilistic Approach

It might be striking at the first glance that econometrics treats all variables and the entire sample as random. 

Let's think about $Y_{i}$ as the salary of individual $i$  (or you can give a name to $i$ if it's too abstract, say Alice). The key here is that before we ask $i$, her salary to **us** is random. We can think about her salary as drawn from some **distribution**. For instance, if we know $i$ is a 22 years old female college graduate coming from China, we may think $Y_{i}$ is some random number drawn from the income distribution of all 22 years old females in China. 

Then what about the numbers in the data set we just saw in [Section 2.2](#2.2 Data Set as A Matrix)? They are real numbers. And the matrix, or data set we saw, is fixed. Where is the randomness? 

Those numbers you saw are **realizations** of random variables $Y,X_{1}$ and $X_{2}$. For each of them, we independently draw $n$ times from their distributions, and observe all the realizations. So we first have a random matrix as follows:
$$
\begin{pmatrix}
Y_{1}&X_{11}&X_{21}\\
Y_{2}&X_{12}&X_{22}\\
Y_{3}&X_{13}&X_{23}\\
Y_{4}&X_{14}&X_{24}\\
\end{pmatrix}.\tag{1}\label{eq.matrix}
$$


All entires of them are random. Then **computer** draws the realizations and put numbers into corresponding entries.

For real world problems where the numbers are not drawn by a computer, we can think about the process in the following way (kind of like a dynamic incomplete information game):

1. Nature (some God) moves first by rolling a dice. This determines the distributions of all the random variables we care about. This dice is in the hand of God. Researchers (you and I) want to make a strategy to back out useful information about the dice. 
   + Income distribution, gender distribution. Mean and variance may be interesting.
   + A bag of red and black balls. Proportion of red balls may be interesting.
2. We first make a plan to collect multiple draws from the distribution. Each draw is an observation (one row in the random sample matrix). All draws of all variables are random.
   + We plan to survey 1000 people asking their income and gender. They are random before we see the results and are from the unknown income distribution set by God.
   + We plan to draw a ball from the bag with replacement for 10 times. The color of the 10 balls drawn from the bag is random before we look, and it follows the unknown distribution set by God.
3. So far, everything is random. With the to-be-seen random draws at hand, we now propose a *plan* to make a good guess of the dice in Step 1. Depending on the purpose, this plan is sometimes called an *estimator*, sometimes a *test*. It is a formula based on the random draws telling us how to calculate the parameters about the dice once we observe the realizations of the random draws.
   + Therefore, **an estimator or a test statistic is also random** because they are functions of random variables.
   + (If you're familiar with the incomplete information game) Here we hope our plan **maximizes the expected payoff**, that is, the random variables may have all possible realizations but before we observe them (i.e., *ex ante*), we want to minimize the cost that our plan leads to a bad guess of God's dice. (If you are familiar with OLS, you can now rationalize why it is defined as the the minimizer of the mean squared error (MSE); MSE is one version of ex ante cost that a researcher faces. )
   + Don't worry if this looks too abstract and does not make much sense. We'll come back to this point again and again when we study concrete estimators.
4. Now we observe the outcomes of the random draws. Immediately, the probability cloud collapses and we see the realizations. 
   + We collect the surveys and see responses of the participants.
   + We open our hands and look at the color of the ball.
5. Regardless of the realizations, we stick to our plan in Step 3. That is, we substitute the realizations into the formula in Step 3. We then get a number which is an **estimate** of the God's dice.

Note that econometric methods we will study kick in as early as Step 3. That is, our plan is not affected by the realizations and it works no matter whether $Y_{1}=100$ or $1000$, whether $X_{1}$ is generated by computer or by survey or by experiment, or whether $X_{2}$ refers to age or GDP. Our methods always work as long as certain assumptions about the dice of God are satisfied.

* Technically, there is a subtle difference between the terms **sample** and **data set**. Sample refers to the random matrix whereas data set refers to the matrix of realized values. Yet in conversations, not everyone follows or pays attention to this difference and we don't have to be this precise in this course as well.

Throughout the semester, we will use and develop tools in probability and mathematical statistics to study econometrics. We will review them after matrix algebra.

### 2.4 Types of Data

From different angles, we can categorize data sets differently.

* Observational data *vs* experimental data. 
  + Observational data contain information that is already there (or, the random variables have already been realized). Someone collects such pre-existing information to make a data set. Survey data (income, education, age, sex), financial data (stock returns), macro data (GDP, CPI) are almost all observational data. Indeed, most economic data are observational.
  + Experimental data are obtained by running experiment. Experiment is more and more popular because causality, the holy grail in economic research, is much easier to be recovered in such data. 
  + We'll study techniques to handle both types of data in the semester.
* Cross-sectional *vs* time series *vs* panel data.
  + Cross-sectional: multiple entities at one time period.
    + A one-time survey of 1000 college graduate. ($n=1000$).
    + Macro performance of G20 in year 2021. ($n=20$).
    + Prices of the 500 stocks in Standard & Poor yesterday. $(n=500)$.
  + Time series: one entity at multiple time periods.
    + Annual GDP of China in the past decade. ($T=10$).
    + Daily price of a specific stock in the past year. ($T=365$).
  + Panel: multiple entities at multiple time periods.
    + A tracking survey of 1000 households for 20 years. 
    + Annual macro performance of G20 in the past decade.
  + We mainly focus on cross-sectional data in this semester, and a bit panel if time allows.

## 3. Causal Inference and Prediction

In some sense, causality and predictions are two primary goals of all braches of data science. Econometrics in particular focuses on causality instead of predictions. This is the most distinct difference between econometrics and stats/ML.

Examples of causal inference:

* Does global warming cause economic recession? (Environment econ)
* Does smoking during pregnancy cause unhealthy babies? (Health econ)
* Does graduate school cause higher income? (Education/labor econ)
* Does minimum wage cause decreased labor demand? (Labor econ)

Example of predictions:

* Forecast tomorrow's stock price based on the price history in the past 3 months. (Finance)
* Forecast next quarter's GDP based on key macro indicators. (Macro)
* TikTok/Bilibili/Amazon/TaoBao push personalized videos/advertisement based on your browsing history. (ML)

As said, econometrics today mainly focuses on causal inference. Like Medical Science, perhaps the most reliable approach to causal inference is to conduct experiments. However, it is not always feasible/legal/moral to conduct experiments in economics (think about the examples above). Fortunately, econometricians have developed a rich toolkit that traces out causality from observational data.

In this semester, we will study both tools based on experimental data and tools based on observational data.
