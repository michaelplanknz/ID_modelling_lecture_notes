
## Motivation and scope

These lecture notes cover the fundamentals of mathematical models of infectious disease dynamics, with a focus on directly transmitted infections in human populations. The course covers the underpinning concepts and basic elements of a range of modelling approaches, with a focus on how different types of model and aspects of the subject are interrelated. We provide references to the literature for further reading, rather than going into great depth in any specialised area. The topics selected are influenced by the author's own research interests and, beyond covering the fundamentals of compartment-based models and branching process models in homogeneous populations, do not attempt to be completely comprehensive.

The course contains a mixture of theoretical results and computational approaches. We approach the subject from a mathematical perspective, but we adopt an expository style that prioritises intuition over rigour, something that will suit some readers and disappoint others. We cover theoretical results primarily where they build conceptual understanding, provide general intuition or insight into infectious disease dynamics, or apply across a range of different model types, rather than purely for their own sake. In other situations, we are content with a computational approach. 

We emphasise interpretation of modelling results and caution that models are just that, and should not be taken too literally. Except for specific cases, such as near-term forecasting, models generally do not are not intended to provide quantitative predictions. Their value lies more in building conceptual understanding of how fundamental mechanisms shape epidemic dynamics at the population scale, comparing the consequences of alternative assumptions and scenarios, uncovering the effects of potential interventions (which can be counterintuitive), and providing a structured framework for the interpretation of epidemiological data. 

As with learning any new subject, doing is the key to understanding. Readers to encourage to verify the algebraic derivations and solutions that are covered, and to experiment with the code blocks provided to explore the effects of changing model parameters and compare the behaviour of different models. 




## Pre-requisites

This course assumes the reader has previously studied the following topics:

1. Calculus of functions of a single variable. 
1. Fundamentals of linear algebra.
1. Ordinary differential equations, including phase plane and linear stability analysis and numerical solutions.
1. Fundamentals of probability, including conditional probability and independence, random variables, and discrete and continuous probability distributions.
1. Implementation of numerical methods and stochastic simulation methods in a computer programming language such as R, Matlab or Python.

In addition, experience with likelihood functions and their use in frequentist or Bayesian inference will be helpful for Topic 8, while familiarity with the behaviour of partial differential equations will be helpful for Topic 11.



## List of topics

Topics 2--7 cover the foundations, while Topics 8 and 9 cover ideas that will underpin any application of models to real-world data. Topics 10--12 cover more advanced or specialised subjects and can be taken selectively. 

1. Introduction. 
1. SIR model.
1. Modelling vaccination.
1. Kermack-McKendrick model.
1. SIRS model.
1. Branching processes.
1. Stochastic models.
1. Model fitting and parameter inference.
1. Reproduction number estimation.
1. Multi-type models.
1. Age-structured models.
1. Network models.

