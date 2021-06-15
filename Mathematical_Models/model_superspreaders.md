\documentclass{article}
\usepackage[utf8]{inputenc}

\title{Model of Super Spreading Events}
\author{Hannah Craddock }
\date{June 2021}

\documentclass[20pt]{article}

\setcounter{secnumdepth}{0}
\usepackage{amsmath,amssymb,amsfonts, tabto}
\usepackage{eucal}
\usepackage{geometry}
%\usepackage[margin=0.5in]{geometry}
%\setlength{\parindent}{0pt}.

 \geometry{
 a4paper,
 %total={160mm,255mm},
 top=19mm,
 bottom = 19mm,
 left = 24mm,
 right = 24mm,
 }


\begin{document}

\maketitle

\section{Model}

A Poisson Mixture model with unknown rates $\lambda_1$ and $\lambda_2$ and mixture weight $p$ is proposed to describe Super-spreading events in an epidemic. 
\\
\\
Likelihood
\\
\\
 $L(\lambda_1, \lambda_2, p | \; \bold{x}) = \displaystyle\prod_{t=1}^{N days} p \cdot [ exp(- \lambda_1_{t})\dfrac{1}{x_t!} \lambda_1_t^{x_t} ]   +  (1-p) \cdot [ exp(- \lambda_2_t) \dfrac{1}{x_t!} \lambda_2_t^{x_t} ] $
\\
\\
where
\\
\\
$\lambda_1_t = r0_1\displaystyle\sum_{i=1}^{t-1}x_i \; ( \; Gamma((t-i); k_1, \theta_1)- Gamma((t-i-1); k_1, \theta_1) $
\\
\\
$\lambda_2_t = r0_2\displaystyle\sum_{i=1}^{t-1}x_i \; ( \; Gamma((t-i); k_2, \theta_2)- Gamma((t-i-1); k_2, \theta_2) $
\\
\\
\\
Log Likelihood
\\
\\
 $l(\lambda_1, \lambda_2, p | \; \bold{x}) = \displaystyle\prod_{t=1}^{N days} log \left( p \cdot [ exp(- \lambda_1_{t})\dfrac{1}{x_t!} \lambda_1_t^{x_t} ] \right)  +  log\left( (1-p) \cdot [ exp(- \lambda_2_t) \dfrac{1}{x_t!} \lambda_2_t^{x_t} ] \right) $
\\
\\
 $l(\lambda_1, \lambda_2, p | \; \bold{x}) = \displaystyle\prod_{t=1}^{N days} log(p) + log(exp(- \lambda_1_{t})) + log(\dfrac{1}{x_t!}) + log( \lambda_1_t^{x_t})  +  log(1-p) + log(exp(- \lambda_2_t))  +log(\dfrac{1}{x_t!}) + log( \lambda_2_t^{x_t}) $
\\
\\
\\
 $l(\lambda_1, \lambda_2, p | \; \bold{x}) = \displaystyle\prod_{t=1}^{N days} log(p) -\lambda_1_{t} + log(\dfrac{1}{x_t!}) + log( \lambda_1_t^{x_t})  +  log(1-p) - \lambda_2_t  +log(\dfrac{1}{x_t!}) + log( \lambda_2_t^{x_t}) $
\\
\\
\\
\textbf{Terms that do not contain one of the model parameters are removed giving};
\\
\\
 $l(\lambda_1, \lambda_2, p | \; \bold{x}) = \displaystyle\prod_{t=1}^{N days} log(p) -\lambda_1_{t} + log( \lambda_1_t^{x_t})  +  log(1-p) - \lambda_2_t + log( \lambda_2_t^{x_t}) $
\\
\\

\subsubsection{Likelihood II}
\\
\\
 $L(r0_1, k_1, \theta_1, r0_2, k_2, \theta_2, p | \; \bold{x}) = \displaystyle\prod_{t=1}^{N days} p \cdot [ exp(- \lambda_1_{t})\dfrac{1}{x_t!} \lambda_1_t^{x_t} ]   +  (1-p) \cdot [ exp(- \lambda_2_t) \dfrac{1}{x_t!} \lambda_2_t^{x_t} ] $
\\
\\
\section{Overview}

Bayesian Epidemic Modelling

\section{Mathematics}

 
% \usepackage[tmargin=1in,bmargin=1in,lmargin=1.25in,rmargin=1.25in]{geometry}.
 

Likelihood
\\
\\
 $L(r0, k, \theta | \; \bold{x}) = \displaystyle\prod_{t=1}^{N days} exp(- \lambda_t)\dfrac{1}{x_t!} \lambda_t^{x_t} $
\\
\\
where
\\
\\
$\lambda_t = r0\displaystyle\sum_{i=1}^{t-1}x_i \; ( \; Gamma((t-i); k, \theta)- Gamma((t-i-1); k, \theta) $
\\
\\
\\
Log likelihood
\\
\\
$l = \displaystyle\sum_{t=1}^{N days} x_t ln(\lambda_t) - \lambda_t$
\\
\\
\\
$l = \displaystyle\sum_{t=1}^{N days} x_t ln(r0\lambda_t) - r0\lambda_t$
\\
\\
\\
$l = \displaystyle\sum_{t=1}^{N days} \left( x_t ln(r0) +  x_t ln(\lambda_t) - r0\lambda_t )\right$
\\
\\
where
\\
\\
$\lambda_t = \displaystyle\sum_{i=1}^{t-1}x_i \; ( \; Gamma((t-i); k, \theta)- Gamma((t-i-1); k, \theta) $
\\
\\
\\
$\deltal l = \displaystyle\sum_{t=1}^{N days} \left( x_t ln(r0\lambda_t) - r0\lambda_t )\right$
\\
\\

\subsubsection{MLE}

$ \dfrac{d l}{d r0} = \displaystyle\sum_{t=1}^{N days}  \dfrac{x_t}{r0} - \displaystyle\sum_{t=1}^{N days} \lambda_t $
\\
\\
\\
$ \hat{r0} =  \dfrac{\displaystyle\sum_{t=1}^{N days} x_t}{\displaystyle\sum_{t=1}^{N days} \lambda_t}$

\subsubsection{Bayesian Inference}
Likelihood
\\
\\
 $L(r0, k, \theta | \; \bold{x}) = \displaystyle\prod_{t=1}^{N days} exp(- \lambda_t)\dfrac{1}{x_t!} \lambda_t^{x_t} $
\\
\\
where
\\
\\
$\lambda_t = r0\displaystyle\sum_{i=1}^{t-1}x_i \; ( \; Gamma((t-i); k, \theta)- Gamma((t-i-1); k, \theta) $
\\
\\
\\
Prior
\\
\\
$p(r0) = Gamma(\alpha, \beta) = \dfrac{\beta^{\alpha}}{\Gamma(\alpha)} r0^{\alpha - 1}e^{- \beta r0} $


\subsubsection{Bayesian Inference}
- Shape, scale** To write out 
\\
Posterior $p(r0|x)$
\\
\\
Given Gamma(1, 1) prior on r0;
\\
\\
 $p(r0|x) \propto  \displaystyle\prod_{t=1}^{N days} exp(- r0\lambda_t)\dfrac{1}{x_t!} (r0\lambda_t)^{x_t} \times r0^{\alpha - 1}e^{- \beta r0} $
\\
\\
\\
$\propto exp(- r0(\displaystyle\sum_{t=1}^{N days}  \lambda_t + 1)) \times (r0\lambda_t)^\(\sum_{t=1}^{N days} x_t $
\\
\\
$\therefore p(r0|x) \propto Gamma(\displaystyle\sum_{t=1}^{N days} x_t + 1, \displaystyle\sum_{t=1}^{N days} \lambda_t + 1)$
\\
\\
$\propto Gamma(\displaystyle\sum_{t=1}^{N days} x_t + \alpha, \displaystyle\sum_{t=1}^{N days} \lambda_t + \beta)$
\\
\\
\\
\textbf{where} (*Above is shape, rate. Below is shape, scale)
\\
\\
$\lambda_t = \displaystyle\sum_{i=1}^{t-1}x_i \; ( \; Gamma((t-i); k, \theta)- Gamma((t-i-1); k, \theta) $
\\
\\
\\



\subsubsection{Bayesian Inference v0}
\\
Posterior $p(r0|x)$
\\
\\
Given Gamma(1, 1) prior on r0;
\\
\\
 $p(r0|x) \propto  \displaystyle\prod_{t=1}^{N days} exp(- r0\lambda_t)\dfrac{1}{x_t!} (r0\lambda_t)^{x_t} \times r0^{\alpha - 1}e^{- \beta r0} $
\\
\\
\\
$\propto exp(- r0(\displaystyle\sum_{t=1}^{N days}  \lambda_t + 1)) \times (r0\lambda_t)^\(\sum_{t=1}^{N days} x_t $
\\
\\
$\therefore p(r0|x) \propto Gamma(\displaystyle\sum_{t=1}^{N days} x_t + 1, \displaystyle\sum_{t=1}^{N days} \lambda_t + 1)$
\\
\\
$\propto Gamma(\displaystyle\sum_{t=1}^{N days} x_t + \alpha, \displaystyle\sum_{t=1}^{N days} \lambda_t + \beta)$
\\
\\
\\
\textbf{where} (*Above is shape, rate. Below is shape, scale)
\\
\\
$\lambda_t = \displaystyle\sum_{i=1}^{t-1}x_i \; ( \; Gamma((t-i); k, \theta)- Gamma((t-i-1); k, \theta) $
\\
\\
\\
Prior
\\
\\
$p(r0) = Gamma(\alpha, \beta) = \dfrac{\beta^{\alpha}}{\Gamma(\alpha)} r0^{\alpha - 1}e^{- \beta r0} $


\section{Explanation of Mathematics}

- Simulation of Epidemic
- Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' = Sum of all people
- Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
-Assumption: Number of daily cases follows a Poisson distribution. Reason why there is spikes from day to day
  
- Metropolis Hastings step
- The logarithm of the acceptance probability $\lambda$ includes; the sum of the log likelihood, prior and the proposal for the current time Y and the previous time step r0[t-1]

- Symmetrical proposal distribution so the proposal for the current time Y and the previous time step r0[t-1] cancel out

\section{Assumptions}
- Number of daily cases follows a Poisson distribution
- Poisson-Gamma conjugacy leads to nice posterior - although not realistic

*
Assume the following Poisson model of two regimes for
n random variables

\end{document}

