
## Introduction to NUNC

NUNC is an acronym for Non-Parametric Unbounded Changepoint detection and is designed to perform the detection of a changes in the distribution of data that has been drawn from an unknown distribution. Furthermore, the algorithm runs in a sliding window. This allows for sequential changepoint detection to take place, along with affording NUNC the ability to deal with non-stationarity in the stream of data. The primary purpose of NUNC is to detect changes in non-stationary or complex distributions of data, where existing parametric detection techniques may fail. The second purpose of NUNC is to perform the analysis online. This means that the data observed by time $t$ is processed before the next data point arrives, and requires that the memory requirements of the algorithm do not overflow as the algorithm is left running and $t \rightarrow \infty$.

The NUNC package contains functions that implement both offline sequential changepoint detection from a pre-recorded dataset, and online changepoint detection where the input is a data generating process rather than a fixed set of data. Furthermore, the package contains implementations of all three variants of the NUNC algorithm. The NUNC algorithms are R wrappers for C++ code, providing fast changepoint detection performance. Furthermore, the online implementation of NUNC utilises the Boost library, and a circular buffer structure, to allow for a stream of data observed in real time to be monitored for a change in distribution.

So that NUNC is able to identify changes in distribution for data that has been drawn from an unknown distribution, a non-parametric approach is taken. This means that NUNC specifies the distribution based only on the observed data, and this distribution changes as different data is observed, rather than relying on classical assumptions such as Gaussianity. As will be discussed in the mathematical background Section, NUNC compares the empirical CDF for the pre and post change data and, if the difference is significant enough, then a changepoint is identified in the data. The comparison takes place using a cost function, and this cost function is formulated by aggregating the difference between the two CDFs at K different quantiles.

An example implementation of NUNC for an online data generating process is presented below. Further examples for both NUNC online, and NUNC offline, are presented later.

```r
set.seed(10)
onlinedata <- c(rnorm(100, 1, 10), rcauchy(100, 12, .1))
 
 f <- function() {
   out <- onlinedata[i]
   i <<- i + 1
out
 }
i <- 1
NUNC_run <- NUNC(f, w = 25, beta = 8, K = 6, method = "local")
#f is the data generating process, w is the window size, beta is the threshold and K the number of quantiles
NUNC_run$changepoint  #changepoint location
#> [1] 200
NUNC_run$t            #time detected (edge of sliding window)
#> [1] 210
```

 <video controls loop><source src="figure/ani1.webm" /><p>plot of chunk ani1</p></video>

This vignette is divided into three parts. The first explores the mathematical background for NUNC, the second presents the NUNC package and the functions it contains, and the final section uses NUNC to analyse a dataset drawn from an accelerometer.

## Mathematical Background for NUNC

NUNC is based upon the comparison of the empirical CDF for independent data points $x_1, ..., x_t$ drawn from a distribution $F_{1:t}$. The empirical CDF, $\hat{F}_{1:t}$, provides an estimate for the distribution of the data, and is given by
$$
\hat{F}_{1:t}(q) = \frac{1}{t} \left\{ \sum_{j=1}^{t} \mathbb{I}(x_j < q) + 0.5 \times \mathbb{I}(x_j = q) \right\},
$$
for some quantile $q$. Using the independence assumption, it can be shown that $t\hat{F}(q) \sim \mathrm{Binom}(t, F(q))$, and a likelihood for a set of data points can then be computed. We write the log-likelihood of the segment $x_{\tau_1+1}, \ldots, x_{\tau_2}$ as
$$
\mathcal{L}(x_{\tau_1+1:\tau_2}; q) = (\tau_2 - \tau_1) \left[\hat{F}_{\tau_1+1:\tau_2}(q) \log(\hat{F}_{\tau_1+1:\tau_2}(q)) - (1-\hat{F}_{\tau_1+1:\tau_2}(q))\log(1-\hat{F}_{\tau_1+1:\tau_2}(q)) \right].
$$

This likelihood can then be used to define a statistic for the detection of a change at a single quantile of the distribution based on the likelihood ratio test. This statistic is defined as
$$
\max_{1\leq \tau \leq t} \mathcal{L}(x_{1:\tau}; q) + \mathcal{L}(x_{\tau+1:t}; q) - \mathcal{L}(x_{1:t}; q),
$$ 
and then in order to increase power we aggregate this statistic across $K$ quantiles. This gives the following test statistic 
$$
C_K(x_{1:t}) = \max_{1 \leq \tau \leq t} \frac{1}{K} \sum_{k=1}^K 2 \left[ \mathcal{L}(x_{1:\tau}; q_k) + \mathcal{L}(x_{\tau+1:t}; q_k) - \mathcal{L}(x_{1:t}; q_k) \right],
$$
and it is this test statistic that is used by NUNC to determine if a change occurs.

NUNC Local is based on this test statistic, although restricts the search to a point inside the sliding window of size W. Furthermore, a point is only identified as a change if the test statistic takes a significant enough value, and this significance is measured by the exceedance of a threshold denoted by $\beta$. The stopping rule for for a change at time $t - W + 1 \leq \tau \leq t$ is
$$
\max_{t-W+1 \leq \tau \leq t} \sum_{k=1}^K 2 \left[ \mathcal{L}(x_{t-W+1:\tau}; q_k) + \mathcal{L}(x_{\tau+1:t}; q_k) - \mathcal{L}(x_{t-W+1:t}; q_k) \right] \geq K\beta.
$$
Using this test statistic, NUNC Local can be seen as identifying changes in the distribution of the data inside the sliding window. 

One drawback of NUNC Local is that $W$ checks are performed inside each sliding window, despite the test statistic being correlated between nearby points and the power of the test being reduced at the edges of the window. To handle this, a subset of the points in the window can be tested instead. This is implemented in the package using the argument \code{max_checks} in the NUNC functions.

In contrast to NUNC Local, NUNC Global does not search for changes inside the sliding window but instead compares the data in the window to the historic data that has been seen so far. In order to store the historic information in a memory efficient manner we again fix $K$ quantiles and update the longrun CDF, denoted by $z^{(t)}_W(\cdot)$, each time a point leaves the sliding window. This update step is denoted by the recursive equation
$$
z^{(t+1)}_W(q) = \frac{1}{t-W+1} \bigg [ (t-W) z^{(t)}_W(q) + F_{t-W+1:t-W+1}(q). 
\bigg ].
$$

This can then be used to compute the likelihood functions for the data, and a change detection test performed in a similar manner to as defined above.

One question that remains is how many quantiles to use, and what value they should take. In line with Haynes (2017) we propose that $K = \lceil a \log(W) \rceil$, and select the value of the quantiles from the sample of data using the formula 
$$
q_k = \hat{F}^{-1} \left(1 + \left(2W+1\right) \exp \left[ \frac{c}{K}(2k-1) \right] \right)^{-1}.
$$
The advantage of selecting the quantiles using this formula is that is gives more weight to the tails of the distribution than would be offered by using quantiles that are evenly spaced.

By selecting the $K$ quantiles using the sample of data, the need to understand the underlying distribution for the data is avoided. Despite this, in some circumstances NUNC may be used to monitor a process with a baseline from a known distribution. In this case, it may be advantageous to specify the quantiles from this known distribution, rather than estimating them from a sample of $W$ points. The reason for this is that the estimates of the tail quantiles may be poor in situations where $W$ is not large enough. It is this issue that NUNC Semi-Parametric addresses, by allowing the user to specify the number of quantiles to use through use of the \code{quantiles} argument.

Having discussed the mathematics underpinning NUNC, we can now explore how the NUNC package can be used to detect changes in the distribution of data. We first present a series of examples that detail the different functions in the NUNC package. After this, we analyse an accelerometer dataset that is contained in the package.

## Detecting Changepoints Using NUNC

Three different forms of NUNC are provided in the NUNC package: NUNC Local, NUNC Global, and NUNC Semi-Parametric. NUNC Local tests for changes in the distribution of the data inside the sliding window, whereas NUNC Global tests for changes in the distribution of the data by asking whether the observations in the window are drawn from the same distribution as the historic observations. In both these methods, the test uses quantiles of the empirical CDF that have been selected using the observed data. In some cases, however, it is argued that there will be a strong belief as to what the distribution of the data is, and so NUNC Semi-Parametric allows for the user to specify the quantiles that are used.

In this section we will explore the different variants, and present code for implementing the algorithms on several different examples. We begin by examining the offline variants of NUNC. For NUNC Local we consider a change in a set of mixture distributions: 

```r
set.seed(1)
x_pre <- c(rnorm(500), rnorm(200, 1, 4), rnorm(300, 2, 2))
x_pre <- sample(x_pre, 1000)

x_post <- c(rnorm(100), rnorm(200, -3, 1), rnorm(200, 3, 2), rnorm(500, -1, 1))
x_post <- sample(x_pre, 1000)

x <- c(x_pre, x_post)

#we use [1:3] to supress the data output

NUNCoffline(data = x, w = 200, beta = 12, K = 4*log(200), method = "local")[1:3]
#> $t
#> [1] 1176
#> 
#> $changepoint
#> [1] 1064
#> 
#> $method
#> [1] "local"
```
We see that NUNCoffline has returned three outputs: \code{t}, the time the change is detected; \code{changepoint}, the location of the change in the window; and \code{method}, which is the method used. We see that a similar analysis can be obtained in a shorter period of time using the \code{max_checks} argument. This is specified as an additional argument using the \code{params} input: this is a \code{list()} containing the argument \code{max_checks =} and can be written as follows:

```r
A <- Sys.time()
NUNCoffline(data = x, w = 200, beta = 12, K = 4*log(200), method = "local")$changepoint
#> [1] 1064
B <- Sys.time()
NUNCoffline(data = x, w = 200, beta = 12, K = 4*log(200), method = "local", 
            params = list(max_checks = 30))$changepoint
#> [1] 1055
C <- Sys.time()

paste("NUNC Local has a", B-A)
#> [1] "NUNC Local has a 0.998739957809448"
paste("NUNC Local on a subset of points in the window has a ", C-B)
#> [1] "NUNC Local on a subset of points in the window has a  0.154743909835815"
```
We next consider analysing a set of data using NUNC Global. We illustrate this using an example of data that is heavy tailed:

```r
set.seed(2)
x_pre <- rcauchy(1000)
x_post <- rcauchy(1000, 1)
x <- c(x_pre, x_post)

NUNCoffline(data = x, w = 100, beta = 5, K = 4*log(100), method = "global")[1:4]
#> $t
#> [1] 1146
#> 
#> $changepoint
#> [1] 1045
#> 
#> $zVals
#>  [1] 0.989005736 0.003824092 0.007648184 0.024856597 0.051625239 0.078393881 0.156787763 0.236137667 0.331739962 0.464627151 0.638623327 0.735181644
#> [13] 0.847036329 0.904397706 0.935946463 0.953154876 0.960803059 0.975143403
#> 
#> $Qhistory
#>    [1]  -4.579796509  -3.269237018   0.193557060  -2.765435136  -5.755257101  -7.580719757  -9.875498015  -3.776484084   0.192672273   8.877668305
#>   [11]   8.903651502   8.832040280   6.347864911   7.236238296   6.610987325   5.619668675   4.940763030   4.403751019   3.104138646   4.582494007
#>   [21]   4.627326040   4.612429770   4.195301712   3.931847587 -10.161138528 -10.545177262 -10.951120916 -10.575746534 -10.129010280  -9.709576721
#>   [31]  -9.371281590  -8.990543740  -8.652261073  -4.556955495  -3.746404714  -5.901634526  -5.722844672  -5.256291792  -5.300344964  -5.411285716
#>   [41]  -5.263463788  -5.573909287  -5.705129652  -5.803989578  -5.802475654  -5.893449361  -5.563165281   8.739041290   8.430220591   7.967304627
#>   [51]   7.298101390   7.074390207   6.339758120   6.471244542   6.112725534   5.806845744   5.258273665   5.543751190   5.451212111   5.086361009
#>   [61]   4.479598294   5.429129417   5.258265369   4.502843813   5.960768975   6.577592084   6.101710101   5.953191797   3.626233799   3.905557212
#>   [71]   4.137947002   3.818130193   3.698280014   3.365082473   3.303680740   3.265878849   3.109290370   3.716943500   2.949828440   3.195858790
#>   [81]   3.238766308   3.133781190   3.034293735   3.010661490   2.962852236   2.968824970   3.158587725   3.181282676   3.268666510   3.309116674
#>   [91]   3.324041346   3.045265250   2.939742468  -4.173674633  -4.069271604  -5.189502233  -4.928410679  -4.640200501  -7.534876697  -7.376297565
#>  [101]  -9.483982272  -8.735818023  -8.509795015  -7.706302064  -8.128604874  -7.779259073 -12.516656829 -12.212195873 -11.660310203 -11.316854557
#>  [111] -10.655170400 -11.724589593 -12.407189126 -12.157434029 -11.869200705 -11.401169094 -11.094006520  -9.966527418 -10.156563446  -9.928132792
#>  [121]  -9.825807778  -7.890233955  -2.895390536  -2.709106293  -2.680814514  -2.384696344  -0.925568816  -0.805624443  -0.899989291  -0.856681173
#>  [131]  -0.674295924  -0.580034166   0.092344957   0.079427301  -0.006081694  -0.423736389  -0.370894069  -0.315101427  -1.386634980  -1.350056610
#>  [141]  -1.347538616  -1.314616590  -1.410171068  -3.360811449  -3.278751883  -2.882673447  -1.772003152  -3.211924030  -2.740318345  -2.624608881
#>  [151]  -3.564064971  -3.016607086  -4.156364908  -3.891059365  -4.266170665  -4.064702248  -3.553646879  -3.309566892  -3.372501357  -2.975984411
#>  [161]  -2.735238750  -1.937724021  -1.126934214  -1.449438906  -1.781532014  -1.529825990  -1.058795684  -0.783482241  -1.750666593  -1.526397882
#>  [171]  -1.274338730  -1.622243233  -0.915719769  -1.943222520  -1.644084654  -2.876230862  -2.148295078  -2.713757342  -4.013595915  -3.566286451
#>  [181]  -3.532722507  -3.435042099  -3.222709550  -3.054951563  -2.994387805  -2.082446170  -1.851486511  -1.676743923  -0.894181749  -0.644271859
#>  [191]  -2.077431152  -1.659500145  -4.862057877  -6.235220666  -5.437768441  -5.416910051  -5.424161925  -5.455428083  -5.254157006  -5.068581384
#>  [201]  -9.314742389  -9.692538123  -9.723194837  -9.808208411  -9.974906639  -9.990292847  -9.701852160  -9.746648797  -9.854919145  -9.834549082
#>  [211] -12.584698634 -12.695184596 -13.171751989 -12.993894699 -13.054025125 -12.969269129 -13.867581896 -13.623459654 -13.681697254 -13.357076450
#>  [221] -13.314478574 -12.989968364  -5.124114066  -5.420006413  -5.529741701  -4.942469046  -2.360630824  -2.960642785  -3.049500706  -3.459497016
#>  [231]  -3.551036414  -3.491512492  -4.426460073  -4.390968443  -4.516828068  -4.616265573  -4.724930246  -4.802797144  -4.427484435  -4.363446318
#>  [241]  -4.575299049  -4.583345065  -4.870381493  -4.297089538  -3.744966725  -3.646535433  -4.323401096  -4.329836261  -4.357390103  -4.355300432
#>  [251]  -4.030719328  -3.981875813  -3.962245628  -4.147239149  -4.044685885  -3.911171585  -3.728729522  -3.688687305  -4.086887811  -4.317387255
#>  [261]  -4.220814079  -4.601812222  -5.035529256  -4.913883332  -4.831759182  -4.829705072  -4.894719468  -4.846801444  -3.859385974  -3.884158401
#>  [271]  -3.871146443  -4.366138626  -4.423446165  -4.342641457   0.198345138   0.583019762   0.254432929   0.598658697   0.174288953  -0.367757917
#>  [281]  -0.362700926  -1.690712559  -1.606113801  -1.717707592  -1.248673977  -0.502464616  -0.969121100  -1.046102076  -1.173735723  -1.416450614
#>  [291]  -1.982629791  -1.815819273  -0.600771606  -0.615054725  -1.507092052  -1.481230275  -1.443142094  -1.415934546   1.004277680   1.257568813
#>  [301]  -1.491492145  -3.419364867  -3.294495972  -3.224829621  -3.099076752  -2.652983118  -2.146268192  -1.859857504  -1.712803725  -0.943585648
#>  [311]  -3.541395757  -1.024889424  -0.902326552  -0.208048954   0.672210850  -0.014815006  -0.710584590  -0.224444287  -0.157340368  -0.989098805
#>  [321]  -1.056040084  -1.070131909  -0.648233897  -1.799227085  -1.926762928  -2.792280952  -3.301580450  -2.500089873  -2.342166726  -2.106774347
#>  [331]  -1.948476728  -1.767927504  -2.480129416  -2.307396305  -1.946973508  -1.771137172  -2.692203560  -2.772301038   1.102951624   2.105857755
#>  [341]   2.250018160   2.310791053   2.917469732   3.100200283   2.159433517   1.040865176   1.810885361   2.220464669  -0.303598918  -0.729760028
#>  [351]  -1.481484376  -0.950815345  -0.817168226  -0.868078418  -0.796073757   0.838964885  -0.352787941  -0.222755694   0.160213023  -0.071528375
#>  [361]  -0.364037251  -0.227377348  -0.289776155  -0.681386484  -1.366876050  -3.274928644  -2.705446185  -3.315142152  -2.977101542  -2.715807852
#>  [371]  -3.001165979   0.019427669  -0.923488525  -1.058103320  -1.664346767  -1.962192263  -1.961224085  -2.030680742  -1.907675562  -5.172616150
#>  [381]   1.195219517   3.666272996   3.265407748   3.324648729  -1.260214990  -1.660530391  -5.624037972  -5.510630727  -5.971402374  -5.965021441
#>  [391]  -5.646790227  -6.005582978  -6.242221203  -5.863774381  -5.636762906  -5.434824463  -5.449440485  -5.314008523  -4.502226726  -5.125112071
#>  [401]  -5.461291639  -5.392859654  -5.149094846  -5.558394757  -5.383439889  -6.356955576 -12.296055032 -12.340836675 -12.196519085 -12.629950581
#>  [411] -13.247495304 -11.099664764 -12.024820927 -12.004444197 -12.232455932 -12.639039081 -13.539571769 -15.073477166 -15.317450059 -15.541857394
#>  [421] -15.605299748 -19.032973668 -19.205970280 -19.716511950 -19.801526124 -19.886261101 -20.740764559 -21.571632392 -21.672911169 -22.080972286
#>  [431] -22.202916742 -22.324516817 -22.306316951 -22.177049093 -22.317039943 -22.129739075 -22.361257582 -30.105819680 -18.255392675 -18.577593343
#>  [441] -18.734122817 -26.200247244 -26.520964758 -27.484149622 -27.881366559 -28.017607429 -27.909861967 -27.814925442 -27.973970387 -28.298330972
#>  [451] -29.294607914 -31.444767579 -32.017612763 -32.178766050 -28.008177303 -27.346790449 -27.327813163 -26.353818591 -26.802174343 -26.917427044
#>  [461] -27.022958716 -26.880162191 -21.903516608 -21.757670521 -22.132912766 -12.775653966 -11.654864848 -11.708556917  -8.954777028  -9.111641222
#>  [471]  -6.531497997  -6.547794303  -6.472754055  -1.830269550  -1.860765072  -4.477663240  -4.487638206  -7.485996421  -7.433308587  -1.939990053
#>  [481]  -2.229456515  -2.015514638  -1.585205651   0.557063259   5.562036583   4.571221636   9.188602039   9.148975194   9.549688213   9.344778426
#>  [491]   7.659869392   5.630554900   5.662352688   5.053499045   5.126721036   4.996457379   5.175775171   5.195741754   4.208213536   3.587486023
#>  [501]   3.103238214   2.614282195   2.724591894   3.752808490   3.595739716   3.860493152   7.152877916   7.184775942   7.250544220   7.205200342
#>  [511]   7.500125758   7.757585307   8.093837617   7.800922829   7.926036227   7.963054371   5.480340126   5.322899511   5.354474488   5.464923179
#>  [521]   5.655008439   5.742629620   5.910283556   5.545489263   5.604127714   5.662724000   5.632595193   5.653559607   5.907479355   6.073658806
#>  [531]   6.599507382   6.648292032   7.053304767   8.190410528   8.722745284  10.384078381  10.462523163  10.474796174  12.757908147  12.872799010
#>  [541]  12.937740106  12.137052579  12.265238329   3.360515268   4.072196650   4.910056518   4.915503049   4.805136420   4.824559666   5.196117522
#>  [551]   5.505231946   1.918159161   0.019474452  -0.009760276   0.894631593   1.011464106   0.931493668   0.859424790   0.872520065   0.940971748
#>  [561]   0.638054874   0.667388961  -0.506053347  -0.424196087  -0.241278776  -0.061491147  -0.206333020  -0.116826753  -0.093456422  -0.172220905
#>  [571]  -0.915105793  -0.843529181  -0.920232949   2.953954068   2.905784615   2.864082064   2.810103849   3.269214935   3.248600023   2.578815375
#>  [581]   2.541005957   2.711811725   2.066904820   3.594748883   3.531867216   4.127429808   2.199982176   2.252022385   2.222915986   2.185605411
#>  [591]   1.814082957   1.300794619   1.287237877   1.334301023   1.058017079   1.098027497   1.618106204   2.224700836   2.366869254   2.575729424
#>  [601]   2.600599183   2.816928311   3.259224960   3.187208479   3.222456672   3.251684082   3.060121310   2.932295036   2.124906641   2.021180647
#>  [611]   2.172100252   1.532203627   1.563146000   1.559510999   1.652763303   1.749359712   1.874327415   3.141374560   3.153798261   2.948926850
#>  [621]   2.818508562   2.872103159   2.897603274   2.773001546   2.878174927   2.895085415   2.820582785   2.857865342   2.872737717   2.797252381
#>  [631]   3.510185699   3.560281428   2.958992405   1.704219392   1.408846327   0.885313340   0.297145370   2.917631679   1.119825496   1.170315032
#>  [641]   1.027071746   1.055997047   1.147312565   1.131600529   1.408556945   1.334911539   1.549581564   3.263480812   3.135949545   3.239199920
#>  [651]   3.673074151   2.832483188   3.123827743   3.111278385   1.531549963   1.792938229   1.959898107   2.551913423   2.553473155   2.447469068
#>  [661]   0.055334122  -0.121143134  -0.138268832   0.803833469   0.648372889   0.497240715   1.306112714   1.422942406   1.375109085   2.062916259
#>  [671]   2.782413086   5.362003792   4.907457235   3.098765491   3.146820545   3.156873318   3.015003027   3.084673659   2.631550834   2.476690939
#>  [681]   1.911546870   1.762229763   2.279176423   2.278844024   2.389720065   2.130552027  -0.618942407  -1.238551249  -0.874529717   0.083851957
#>  [691]  -0.100454015   0.569637151   0.453738602   0.616995061   0.849559159   1.034405392   0.681141372   0.596044689  -1.221567390  -2.147200033
#>  [701]  -1.822294819  -1.447433472  -1.918993108  -1.241581568  -1.303613387  -1.450527963  -0.799021154  -1.001441639  -1.364576744  -2.260125553
#>  [711]  -2.864682047  -3.148703306  -2.808311935  -2.683199605  -3.953908408  -3.821534881  -4.095243411  -2.684811960  -2.238151412  -0.699442498
#>  [721]  -0.707870197  -0.706584556  -0.267493739   0.689889056   0.841150709   1.182809073   1.168398205  -0.505755337  -0.318926352  -0.338051767
#>  [731]  -0.456777575  -0.207699521  -0.125086300  -0.208756814   2.755174113   3.445583213  -0.783112884  -2.664946877  -2.387598024  -1.294445735
#>  [741]  -0.200718650  -0.195031670  -0.189345862  -0.168383723   0.945562481   0.933705823   0.504447940   0.203307881  -0.145507756  -0.400916984
#>  [751]  -0.525031880  -0.691623981  -0.556390306  -0.481833681  -0.912543463  -1.714673281  -1.670166634  -1.286951070  -1.044395966  -1.495487832
#>  [761]  -1.036628786  -1.038539029  -1.012397352  -1.537559117  -1.775623722  -2.330148767  -2.287659119  -0.586505269  -0.194012183  -0.460718151
#>  [771]  -0.448393737  -1.718606297  -2.087003013  -2.793162251  -2.660812625  -0.206530448   1.521901184   1.986702924   2.598578440   1.795505754
#>  [781]   1.037724444   1.784324602   1.513238367   1.532907981   1.628546271   1.925068601   2.328593244   2.468406761   2.403863761   2.283968740
#>  [791]   2.332643974   2.732666734   2.649025840   2.680757124   2.561683390   2.093324733   2.206583429   2.203691481   1.762088690   2.095938468
#>  [801]   0.793941491   2.909854285   2.650277083   2.354512430   2.246478340   2.246253553   1.974195793   2.102772720   1.909636544   2.445896296
#>  [811]   2.588286777   2.562074150   2.732051247   2.737711192   4.522899512   4.462923108   3.980551554   3.975891322   4.238228417   1.517368320
#>  [821]   1.380047864   1.395213388   1.815635191   0.959446763   0.991506622   1.114024869   1.156618219   2.347564471   2.817962855   2.910225405
#>  [831]   2.913985072   3.074839410   2.427393041  -0.097738477  -0.876132795  -0.894698592  -0.913230582  -1.446053398  -1.518200050  -1.517325890
#>  [841]  -1.466449988  -1.494777497  -1.534247550  -1.614599974  -1.665092206  -1.671449999  -1.700441251  -1.433994992  -0.881085373  -0.874829029
#>  [851]  -0.486965163  -0.870173362  -0.925779235  -0.971235425  -1.007205908  -1.033351091  -6.616173204  -7.213334871  -3.025453856  -3.061601560
#>  [861]   0.722480726   1.371641244   1.865241127   1.856382842   1.871889694  -0.554453247  -0.570841600  -1.251167492  -0.698537622   1.691062664
#>  [871]   2.132406324   2.031303114   2.020745988   2.497158906   2.036527690   4.373567771   3.369384524   3.261037262   3.208075982   3.199794449
#>  [881]   3.304135597   2.992488745   3.388673473   3.405204475   3.050675143   2.546358710   2.200451567   2.169293078   0.689482692  -1.493160594
#>  [891]  -1.512706315  -1.655253166  -1.727094395  -1.719891166  -1.737388218  -2.481049383  -2.514022042  -2.367276050  -0.362213932  -0.624130374
#>  [901]   0.259679782   1.168926384  -3.669696649  -3.346774329  -3.285533894  -3.716760524  -3.810138741  -3.765640423  -3.743755079  -3.349586055
#>  [911]  -3.668400746  -3.695267514  -4.242763729  -4.368133621  -4.742086430  -4.772149524  -5.138465958  -5.113214383  -5.085696928  -1.638055656
#>  [921]  -1.409212326  -1.101728457  -1.153034549  -0.758604355  -0.654661823  -0.660110732  -0.666190664  -1.067483895  -0.176716777   0.273258781
#>  [931]  -0.473218862   0.105685181   0.470464482   2.323373418   4.842325400   0.132336546  -0.263591277  -0.358396936   1.014592313  -1.031512159
#>  [941]  -1.318589194  -1.216346502  -1.311917210  -1.366836873  -4.067117751  -5.328058488  -5.568920078  -6.248655896  -5.729208907  -5.627507854
#>  [951]  -5.270616776  -4.378151415  -3.822929780  -4.950444234  -5.331464823  -5.570251248  -1.002955279  -0.724679743  -1.305072616  -1.898393639
#>  [961]  -2.275888143  -2.891052088  -2.149266446  -2.197104168  -3.766869365  -0.974053709  -1.050915966  -3.085447208  -3.914683147  -5.085609050
#>  [971]  -4.233885502  -4.404102192  -3.985742828  -2.880936201  -0.967727860  -1.028295794  -1.111367399  -1.166751656   1.129551690  -0.923740167
#>  [981]  -7.313341513  -8.807162913  -7.820005490  -7.728659143  -5.005042871  -5.105694501  -6.542050140  -3.816395688  -0.059555725   5.594283441
#>  [991]   5.806934539   4.720345917   5.462632750   4.959288378   3.852057095   4.894340411   4.474813975   7.859652553   4.753177249   3.148814189
#>  [ reached getOption("max.print") -- omitted 45 entries ]
```
We see that the input structure is similar, however the output structure also contains an argument called \code{zVals}; this is the value of the longrun CDF at the $K$ quantiles for the data.

The third form of NUNC is NUNC Semi-Parametric. We explore how this can provide better results than NUNC Local in cases where the quantiles for the data are pre-specified correctly by examining a scenario where there is a high variance in the distribution for the data under the null hypothesis. This scenario makes estimation of the quantiles from a small sample of data difficult.

As in the \code{max_check} case, the additional \code{quantiles} argument is specified using the \code{params} argument. We see in the following example we obtain a correct changepoint when using NUNC Semi-Parametric, whereas an incorrect changepoint is returned using NUNC Local.

```r
set.seed(3)
x_pre <- rnorm(1000, 0, 10)
x_post <- rnorm(1000, 10, 4)
x <- c(x_pre, x_post)

NUNCoffline(data = x, w = 50, beta = 6, K = 4*log(50), method = "local")[1:3]
#> $t
#> [1] 97
#> 
#> $changepoint
#> [1] 58
#> 
#> $method
#> [1] "local"

quantiles <- qnorm(seq(0.05, 0.99, length = ceiling(4*log(50))), 0, 10)

NUNCoffline(data = x, w = 50, beta = 6, K = 4*log(50), method = "semi-param",
            params = list(quantiles = quantiles))[1:4]
#> $t
#> [1] 1038
#> 
#> $changepoint
#> [1] 987
#> 
#> $zVals
#>  [1] 0.04352227 0.10627530 0.17611336 0.23886640 0.30364372 0.37246964 0.43016194 0.47773279 0.53846154 0.61234818 0.67408907 0.74089069 0.80060729
#> [14] 0.86842105 0.92408907 0.98987854
#> 
#> $method
#> [1] "semi-param"
```
Online analysis of a data stream can also be performed using NUNC, and this is illustrated by the following examples using the \code{NUNC} function. This takes a similar set of inputs to the \code{NUNCoffline} function, however rather than take a \code{data =} argument, we instead use {f =}, where \code{f} is a data generating function that returns a single number at each time point.


```r
 set.seed(4)
 onlinedata <- c(rnorm(300, 1, 10), rcauchy(1200, 12, .1))
 
 f <- function() {
   out <- onlinedata[i]
   i <<- i + 1
    #uncomment and run with x11() or windows() to see real time plotting
    #plot(tail(onlinedata[1:i], 80), ylab = "test data", type = "l")
out
 }
i <- 1
NUNC(f, 50, beta = 10, K = 15, params = list(max_checks = 3), method = "local")
#> $t
#> [1] 324
#> 
#> $changepoint
#> [1] 299
#> 
#> $method
#> [1] "local"
i <- 1
NUNC(f, 50, beta = 10, K = 15)
#> $t
#> [1] 345
#> 
#> $changepoint
#> [1] 294
#> 
#> $zVals
#>  [1] 0.97457627 0.03050847 0.03728814 0.05423729 0.11864407 0.22033898 0.31864407 0.43389831 0.58305085 0.75593220 0.85423729 0.90847458 0.92881356
#> [14] 0.95254237 0.96949153
#> 
#> $method
#> [1] "global"
```
In the above example, the data is fed into the function one point at a time and analysis performed. If the code is run with the plotting aspect of the function uncommented, a user will also be able to see the stream of data plotted whilst it is being observed.

## Analysing Accelerometer Data Using NUNC

Now that a user is familiar with the different functions contained in the NUNC package, we can explore a set of data created using an accelerometer. This dataset is contained within the package, and contains 100 instances where the movement of a Dualshock4 controller has been recorded. The movement is captured by a change in the y axis value recorded by the motion sensor in the controller. 


```r
data(accelerometer)
```

Each entry in the list is a data frame containing three variables:
\describe{
 \item{y}{A vector of length 2000 indicating the y axis position recorded by the controller.}
 \item{changepoint}{The time the action occurs.}
 \item{action}{The action that is performed.}
}
The four actions that can take place are: picking up the controller, sitting with the controller and moving it, sliding the controller along a surface, and shaking the controller. In every example, the movement of the controller changes the distribution of y axis values.

We can explore an example of this change in distribution:
![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

We can now use NUNC to detect changes in the distribution of the data. For example


```r

accelerometer[[52]]$changepoint
#> [1] 497

NUNCoffline(data = accelerometer[[52]]$y, w = 100, beta = 17, K = 10, method = "local",
            params = list(max_checks = 30))$changepoint
#> [1] 573


NUNCoffline(data = accelerometer[[52]]$y, w = 100, beta = 17, K = 10, method= "global")$changepoint
#> [1] 535
```
In both instances of NUNC, there is some delay between the emergence and the detection of the change. This is to be expected given the sequential nature of the algorithm. We also observe the better performance of NUNC Global, which again is to be expected given the fact that the pre-change distribution is stationary. 

We are then able to plot the output of this analysis using the \code{NUNCplot} function, with the changepoint location indicated by the red line:

```
#> Error in NUNCplot(local_run, xlab = "Time", ylab = "Data", main = "NUNC Local"): Must input an object of class NUNCout
#> Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
#> Error in NUNCplot(global_run, xlab = "Time", ylab = "Data", main = "NUNC Global"): Must input an object of class NUNCout
#> Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

We can build upon this and perform a more detailed analysis of the accelerometer dataset using NUNC. In particular, we can assess both the detection power, and the detection delay, for both NUNC Global and NUNC Local. Furthermore, the detection delay is given by for NUNC Local and NUNC Global respectively, along with the standard deviation of this delay. For both variants of NUNC, we use the threshold given by a theoretical result presented in Austin (2021). This is set to control the probability of false alarm under the null hypothesis at the 10% level.

The results for this experiment are given in the table below.


Table: Table showing results from NUNC analysis of accelerometer dataset

|           | Power| Delay| Delay Sd|
|:----------|-----:|-----:|--------:|
|Local      |  0.72| 69.65|   136.38|
|Local Grid |  0.72| 74.53|   142.86|
|Global     |  0.74| 83.31|    94.64|

The code to implement this example is given below:


```r
NUNC_local_detected <- NUNC_grid_detected <- NUNC_global_detected <- 0
NUNC_local_delay <- NUNC_grid_delay <- NUNC_global_delay <- rep(NA, 100)

#compute theoretical thresholds for the test to control false alarms at 10%

local_threshold <- find_beta(alpha = 0.1, t = 2000, w = 100, K = 15, method = "local")
global_threshold <- find_beta(alpha = 0.1, t = 2000, w = 100, K = 15, method = "global")

for(i in 1:100){
  
  local_run <- NUNCoffline(accelerometer[[i]]$y, 100, local_threshold, 15, "local")
  
  grid_run <- NUNCoffline(accelerometer[[i]]$y, 100, local_threshold, 15, "local",
                           params = list(max_checks = 30))
  
  global_run <- NUNCoffline(accelerometer[[i]]$y, 100, local_threshold, 15, "global")
  
  #if the changepoint is detected, and is not before the accelerometer is moved,
  #record a successful detection
  
  if (local_run$changepoint != -1 & local_run$t >= accelerometer[[i]]$changepoint){
    NUNC_local_detected <- NUNC_local_detected + 1
    NUNC_local_delay[i] <- local_run$t - accelerometer[[i]]$changepoint
  } 
  
  if (grid_run$changepoint != -1 & grid_run$t >= accelerometer[[i]]$changepoint){
    NUNC_grid_detected <- NUNC_grid_detected + 1
    NUNC_grid_delay[i] <- grid_run$t - accelerometer[[i]]$changepoint
  } 
  
  if(global_run$changepoint != -1 & global_run$t >= accelerometer[[i]]$changepoint){
    NUNC_global_detected <- NUNC_global_detected + 1
    NUNC_global_delay[i] <- global_run$t - accelerometer[[i]]$changepoint
    }
  
  print(i)
  
}

#we can now compute the power and average detection delay

local_power <- NUNC_local_detected / 100
grid_power <- NUNC_grid_detected / 100
global_power <- NUNC_global_detected / 100

local_delay <- mean(NUNC_local_delay, na.rm = TRUE)
local_delay_sd <- sd(NUNC_local_delay, na.rm = TRUE)

grid_delay <- mean(NUNC_grid_delay, na.rm = TRUE)
grid_delay_sd <- sd(NUNC_grid_delay, na.rm= TRUE)

global_delay <- mean(NUNC_global_delay, na.rm = TRUE)
global_delay_sd <- sd(NUNC_global_delay, na.rm = TRUE)

```
