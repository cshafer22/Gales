ABOUT

Hello and welcome to Gales!

This is the beginning stages of a probability and stochastic processes library
which will grow to include statistics and machine learning capabilities through
OCaml libraries. It utilizes the Zarith library and is built on top of Z and Q
for large int capabilities. At the moment it can accurately compute probabilities
and expectations of many common random variables, both discrete and continuous.
Soon it will grow to encomapss many more.

In order to use please download and install the Zarith library here:
https://github.com/ocaml/Zarith

Once installed, open utop in terminal. To use Zarith library in utop use the commands:

	#require "zarith";;
	#require "zarith.top";;

If there are issues, it is possible to go through topfind to examine possible errors
and listing all available packages:

	#use "topfind";;
	#list;;


CONSTRUCTION OF RANDOM VARIABLES

Consider a normal discrete uniform random variable, we will call "z". We can
build "z" by calling:

	let z = RV.disc_unif 4;;

This will establish the random variable "z" as defined over the float set
S = {1.0, 2.0, 3.0, 4.0}. Notice that if we call the function:

	disc_pr z "x";;

Where x is a member of S the resulting out put will be 0.25, or 1/4. Where
x not in S, the output will be 0.0. We can then take the expectation of z
by calling:

	e z;;

The output will be 2.5. Further we can consider defining a normally distributed
continuous random varible "w" by:

	let w = RV.gauss 0.0 1.0;;

With this we create a Gaussian of mean 0 and standard deviation of 1. Convention
then tells us we can break this bell burve into sections. Indeed consider:

	cts_pr w 0.0 1.0;;
	cts_pr w 1.0 2.0;;
	cts_pr w 2.0 3.0;;
	cts_pr w 3.0 4.0;;

The respective results will be 0.341344746068542926, 0.135905121983277866,
0.0214002339165491051, and 0.00131822678979698349 corresponding to the space
under the bell curve respective to those probability integrals. Of course
summing them together we get a number that approaches 0.5. Further if we
check:

	cts_pr w 0.0 infinity;;

We clearly see that 0.5 matches the output of half the bell curve by design.
Further we can find the expectation of w corresponding to the mean we defined
earlier with:

	e w;;

I encourage all users to play around with the available random variables already
established. Full functionality exists for the current list:

RV.disc_unif (n: int)	    		(discrete uniform random variable)
RV.bern (p: float)			(bernoulli random variable)
RV.binom (n: int) (p: float)		(binomal random variable) 
RV.geo (p: float)			(geometric random variable)
RV.pois (lambda: float)			(poisson random variable)
RV.cts_unif (a: float) (b:float)	(continuous uniform random variable)
RV.exp (lambda: float)			(exponential random variable)
RV.gauss (mu: float) (sig2: float)	(normal/gaussian random variable)

More functionality is on the way! Thank you and enjoy!
