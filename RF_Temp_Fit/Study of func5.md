# Study of $func5$

## Definition

The following function $fun5$ is defined as follows
$$
func5 = A\frac{t-t_0}{(C+|t-t_0|^k)^{(1/k)}} + B
$$
For convenience, in the RF heating problem, it is possible to express also $func5_to10$ as
$$
func5\_to10 = 10^{func5} = 10^{A\frac{t-t_0}{(C+|t-t_0|^k)^{(1/k)}}}10^B
$$
What I do usually is fit the Temperature with $func5\_to10$, and then display both the temperature T and $func5\_to10$ in `semilogy`. The function $func5\_to10$ is plotted below in linear and `semilogy` scales. 

![testplotfunc5_to10](/home/adrien/Documents/Programmes/Python/THESE/RF_Temp_Fit/20220617_data/testplotfunc5_to10.png)

*Additional information can be found [over Wikipedia](https://en.wikipedia.org/wiki/Sigmoid_function). The function $func5$ is a quite general form of a sigmoid function, that can be used to reproduce several well-known functions such as the logistic function, the error function or the hyperbolic tangent function, but it is not the more generalised case because apparently it is possible to define a sigmoid as a special case of another function that can be found in the aforementioned wiki.*



I am now interested on the effect of $B$, $C$ and $k$ parameters. I will present one graph for each parameter where several curves are plotted with all parameters constant, except the chosen one.

## Parametric study

In the aforementioned wiki, one [scientific reference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4692073/) (Dunning et al.) can be found, where a simplified version of the $func5$ function is used and studied. They use a specific case of $func5$ where $B=0.5$ and $C=1$. They call their function $\pi(\beta(t-\alpha))$ with $\beta$ playing the same role as $2A$ and $\alpha$ playing the same role as $t_0$. It is mentioned that

> [k] is a shape parameter governing how fast the curve approaches the asymptotes for a given slope at the inflection point. When *κ*=1 the function is the absolute sigmoid function and when *κ*=2 the function is the square root sigmoid function; when *κ*=1.5 the function approximates the arctangent function, when *κ*=2.9 it approximates the logistic function, and when *κ*=3.4 it approximates the error function. *κ* is estimated in the likelihood maximization. Models with generalized symmetrical protection curves were fitted to each of the illustrative datasets.

In Dunning et al., this $\pi$ function is symmetrical. To generate non-symmetrical functions they refer to a $\Phi$ function and use it in a similar manner referred in the Definition section (see italics remark about the more general case).

### Preliminary remarks

It turns out our $func5$ function can’t produce non-symmetrical curves, which is a problem considering the temperature evolution in the RF heating is not symmetrical in `semilogy`. This is visible when plotting the time derivative of the temperature and $func5$. While the first derivative is non-symmetrical, the derivative of $func5$ remains symmetrical. This problem may be considered in the future if one wants to improve the fit quality and the derivative exactitude.

### Parameter B

The parameter B is the value of the function $func5\_to10$ at $t=t_0$ so that
$$
func5\_to10(t=t_0) = B
$$
![func5_to10_logy_B](/home/adrien/Documents/Programmes/Python/THESE/RF_Temp_Fit/20220617_data/func5_to10_logy_B.png)

### Parameter C

![func5_to10_logy_C](/home/adrien/Documents/Programmes/Python/THESE/RF_Temp_Fit/20220617_data/func5_to10_logy_C.png)

### Parameter k

![func5_to10_logy_k](/home/adrien/Documents/Programmes/Python/THESE/RF_Temp_Fit/20220617_data/func5_to10_logy_k.png)

![func5_to10_logy_k_bis](/home/adrien/Documents/Programmes/Python/THESE/RF_Temp_Fit/20220617_data/func5_to10_logy_k_bis.png)

## Beyond what is already done

The excellent tool to explore Growth-curve can be found [here](https://ogarciav.github.io/grex/). Some cases of non-symmetrical curves could be found with assymetrical derivatives. See Schumacher, or also on the top left part of the clickable graph.