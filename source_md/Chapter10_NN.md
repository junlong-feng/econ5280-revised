# ECON5280 Chapter 10 Neural Networks

<font size="5">Junlong Feng</font>

## Outline

* Motivation: A glimpse of the frontier of deep learning+econometrics.
* Neural Networks and Deep Learning: A powerful way to approximate functions.
* A Global Method to Approximate Conditional Expectations: All the causal parameters in this semester are identified in terms of conditional expectations.
* Applications and Implementation.

## 1. Neural Networks

Neural networks are the building block of deep learning. They have many uses, but here we focus on function approximation, which is important for causal inference: In this semester, we learned 5 causal parameters, ATE, CATE, LATE and CLATE. Two of them, CATE and CLATE, are identified as multidimensional functions. ATE, though a scalar, is identified as a functional of some multidimensional functions. Previously, we used tree-based methods to approximate these functions. Now we introduce neural nets as an alternative approach. 

- Origins: Algorthims that try to mimic the brain.
- Widely used in 80s and early 90s. Less popular in late 90s.
- Today: State-of-art technique for many applications.
- Belongs to the global methods we introduced before.

Example 1. Suppose we have a real vector $x=(x_{1},x_{2},x_{3})'$ and a real valued function of it $h(x)$. We observe the function value $y$, but don't know $h$ is. We want to use a flexible enough **parametric** function $h_{\theta}(x)$ to approximate the unknown $h(x)$. A neural net is such an approximation:



<img src="/Users/junlong/Documents/HKUST Teaching/Fall 2024/NeuralNets.png" alt="Network" style="zoom:25%;" />

<font size="2">Figure from Andrew Ng.</font>

- Each circle is a neuron, or a unit.

- The first layer has 4 units, consisting of the three variables and a constant $x_{0}=1$.

- They are going to be combined linearly (called the input function) and passed to the hidden layer by a known function $f$ (called activation function).

  - $a_{1}^{(2)}=f^{(1)}(\theta_{01}^{(1)}+\theta_{11}^{(1)}x_{1}+\theta_{21}^{(1)}x_{1}+\theta_{31}^{(1)}x_{3})$, $a_{2}^{(2)}=f^{(1)}(\theta_{02}^{(1)}+\theta_{12}^{(1)}x_{1}+\theta_{22}^{(1)}x_{1}+\theta_{32}^{(1)}x_{3})$, etc.
  - Then from the hidden layer to the final outcome, we again linear combine $a^{(2)}$s first by some $\theta^{(2)}$s, and then pass into the second activation function $f^{(2)}$. That is our final $h_{\theta}(x)$:

  $$
  \begin{align*}
  h_{\theta}(x)=&f^{(2)}(\theta^{(2)}_{0}+\theta^{(2)}_{1}a^{(2)}_{1}+\theta^{(2)}_{2}a^{(2)}_{2}+\theta^{(2)}_{3}a^{(2)}_{3}).
  \end{align*}
  $$

  - One example is Let $f(z)=e^{z}/(1+e^{z})$. Called **logistic** function in stats/econometrics and **sigmoid** function in ML.

- So $h_{\theta}$ in this example is a known function up to $\theta$, where $\theta$ contains $(3+1)\times 3+(3+1)\times 1$ parameters.

- Generalize this example:

  - Hidden layer can have more units than the inputs.
    - A single hidden layer neural network with a linear output unit can approximate any continuous function arbitrarily well, given enough hidden units (Hornik (1991), "Approximation capabilities of multilayer feedforward networks", *Neural Networks*).
  - There can be more than one hidden layers.
    - One or two hidden layers: Shallow neural network.
    - Three or more hidden layers: Deep neural network. Deep learning.
  - The outcome can be multiple (pattern recognition etc.) 

Example 2. Sigmoid activation function to approximate the max function. To understand how $h_{\theta}$ can approximate $h$, suppose we have $x_{1}\in \{0,1\}$ and $x_{2}\in \{0,1\}$. Suppose $h(x)=\max\{x_{1}, x_{2}\}$. We can use a sigmoid function $h_{\theta}(x)=\exp(\theta_{0}+\theta_{1}x_{1}+\theta_{2}x_{2})/(1+\exp(\theta_{0}+\theta_{1}x_{1}+\theta_{2}x_{2}))$ to approximate it by setting $\theta_{0}=-20$, $\theta_{2}=40$ and $\theta_{3}=40$. Now we can verify that

- $h_{\theta}(0,0)=\exp(-20)/(1+\exp(-20))\approx 0=\max\{0,0\}$.
- $h_{\theta}(0,1)=h_{\theta}(1,0)=\exp(20)/(1+\exp(20))\approx 1=\max\{0,1\}=\max\{1,0\}$.
- $h_{\theta}(1,1)=\exp(60)/(1+\exp(60))\approx1=\max\{1,1\}$.

Example 3. Play with neural nets at https://playground.tensorflow.org.

## 2. Global Approximation of Conditional Expectation

### 2.1 Conditional Expectation as the Best Predictor

We first need to know the following fact: If we want to find a function of $W$, $f(W)$, to approximate random variable $Y$ in the sense that

$$
\min_{f}\mathbb{E}\left[(Y-f(W))^{2}\right].
$$
Then the minimizer is $f^{*}(W)=\mathbb{E}(Y|W)$.

*Proof.* 
$$
\begin{align*}
\mathbb{E}\left[(Y-f(W))^{2}\right]=&\mathbb{E}\left[(Y-\mathbb{E}(Y|W)+\mathbb{E}(Y|W)-f(W))^{2}\right]\\
=&\mathbb{E}\left[(Y-\mathbb{E}(Y|W)^{2}\right]+\mathbb{E}\left[(\mathbb{E}(Y|W)-f(W))^{2}\right]\\
&+2\mathbb{E}\left[(Y-\mathbb{E}(Y|W))\cdot(\mathbb{E}(Y|W)-f(W))\right]\\
=&\mathbb{E}\left[(Y-\mathbb{E}(Y|W)^{2}\right]+\mathbb{E}\left[(\mathbb{E}(Y|W)-f(W))^{2}\right].
\end{align*}
$$
Hence, the expectation on the left hand side is minimized by setting $f(W)=\mathbb{E}(Y|W)$. 

Therefore, for instance, we are interested in $\mathbb{E}(Y|D=d,W)$, an important component of the identification equation for CATE and ATE, the above result implies that it satisfies:
$$
\mathbb{E}(Y|D=d,W)=\arg\min_{h_{d}}\mathbb{E}\left[\left(Y-h_{d}(W)\right)^{2}|D=d\right].
$$


### 2.2 Putting Everything Together

Although $h_{d}(W)$ is unknown, we can approximate it by neural nets: $h_{d;\theta}(W)$. The minimization problem in (2) thus can be approximated by
$$
\mathbb{E}(Y|D=d,W)\approx\arg\min_{\theta}\mathbb{E}\left[\left(Y-h_{d;\theta}(W)\right)^{2}|D=d\right].
$$
Next, recall that expectation can be approximated by sample average. Hence,
$$
\widehat{\mathbb{E}}(Y|D=d,W)\equiv\arg\min_{\theta}\frac{1}{|\{i:D_{i}=d\}|}\sum_{i\in\{i:D_{i}=d\}}\left(Y_{i}-h_{d;\theta}\left(W_{i}\right)\right)^{2}.
$$
You can estimate, for instance, the propensity score $\Pr(D=1|W)$ in a similar way.

Pros of neural nets:

- Like causal forest, neural nets can handle high dimensional $W$ (that is, a large number of covariates can be included into the model), yielding rich heterogeneity.
- Variables are not necessarily to be numerical.
- GPU computation valid.

Cons:

- Asymptotic inference is still unclear.
  - However, recent research has established that the rate of convergence of the neural net-based estimators is sufficiently fast. So, although the asymptotic distribution of neural nets is unclear, implying that we cannot do inference for CATE estimated by neural nets, we can use neural nets to construct DML-AIPW, and thus valid inference is available for the resulting estimator of ATE.


## 3. Applications and Implementation

In principle, neural nets can be applied to any scenario when you need to estimate conditional expectations. One application is "search-advertisement position effects" in Hartford, Lewis, Leyton-Brown, and Taddy (2017).

- Sponsored advertisements paid to search engine generate over 40 billion U.S. dollars annually for search engines (Goldman and Rao, 2014, "Experiment as instruments: Heterogeneousgeneous position effects in sponsered search auctions", *EAI Endorsed Transactions on Serious Games*). Example: search anything, like health insurance in Google. There can be **Ad** showing up on the first a few results 

- They examine how advertiser's position on the Bing search page (called *slot*, our $D$) affects probability of a user click ($Y$).
- $D$ is endogenous because it is correlated with latent user intent. For instance, when a user searches for "Coke", it is likely that both she will click and that "Coke" is the top-advertiser.
- They use neural nets for nonparametric instrumental variable regression, called *DeepIV*.
- Instruments $Z$ contain a large number of indicators for experiments run by Bing in which advertiser-query pairs were randomly assigned to different algorithms that in turn scored advertisers differently in search auctions, resulting in random variation in position.
- Covariates $W$ contain text data like user query.
- $n>20$ million.
- They find that, among other results, "on-brand query" (like Coke bids on the word "Coke"), position is worth more to small websites.

For implementation, Hartford, Lewis, Leyton-Brown, and Taddy (2017) developed a Python package called DeepIV.

- https://github.com/jhartford/DeepIV.
- The package builds on Keras.
- Keras is available in R.



