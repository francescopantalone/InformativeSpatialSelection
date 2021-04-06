


# Some questions

1. Variogram/Covariogram

* Q1: The formulas in the paper are for the general case <img src="/latex/tex/6bf7dc456eb8f6f7af59197f86dbf30d.svg?invert_in_darkmode&sanitize=true" align=middle width=51.82233374999999pt height=27.91243950000002pt/> or for the case <img src="/latex/tex/ebfa6eaf066249e64c3546ae3b55789b.svg?invert_in_darkmode&sanitize=true" align=middle width=51.531803399999994pt height=26.76175259999998pt/>?

R1: For the definition of <img src="/latex/tex/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92464304999999pt height=22.465723500000017pt/> and <img src="/latex/tex/9b325b9e31e85137d1de765f43c0f8bc.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92464304999999pt height=22.465723500000017pt/>, this is the general case: <img src="/latex/tex/64cd4f64e03c32384b447158216c056a.svg?invert_in_darkmode&sanitize=true" align=middle width=53.648827649999994pt height=27.91243950000002pt/>

2. **Equation (5)**

<img src="/latex/tex/d70551bde9a9ed1fd53007724e859110.svg?invert_in_darkmode&sanitize=true" align=middle width=128.0475042pt height=24.65753399999998pt/> Gaussian random process

<img src="/latex/tex/a41ed0ded6ecb7929d999ec3466144c2.svg?invert_in_darkmode&sanitize=true" align=middle width=79.06582859999999pt height=22.465723500000017pt/>

<img src="/latex/tex/2e27c36555566c8987840e02be9521b5.svg?invert_in_darkmode&sanitize=true" align=middle width=67.60687229999999pt height=22.465723500000017pt/>, <img src="/latex/tex/b78b42acc7ceb812389f4550a1d8ab89.svg?invert_in_darkmode&sanitize=true" align=middle width=118.19109719999999pt height=30.984656999999984pt/>

* Q2: <img src="/latex/tex/b97c06db6ff28ccde87a1b8ff1c3afa9.svg?invert_in_darkmode&sanitize=true" align=middle width=46.21760714999999pt height=22.831056599999986pt/>?

R2: not necessarily.

Expected value of signal <img src="/latex/tex/ed6a9a9a3debf2ed31b3d31ceba8c433.svg?invert_in_darkmode&sanitize=true" align=middle width=74.52783854999998pt height=22.465723500000017pt/>, <img src="/latex/tex/d276a850afef18e55a40ba7e827080d8.svg?invert_in_darkmode&sanitize=true" align=middle width=94.38332474999999pt height=24.65753399999998pt/> and covariance matrix between points <img src="/latex/tex/81948cc655401042d29788c016a4d2c1.svg?invert_in_darkmode&sanitize=true" align=middle width=71.9841309pt height=14.611878600000017pt/> and <img src="/latex/tex/2fd2161c9fff6652d05d0a2dd556da3d.svg?invert_in_darkmode&sanitize=true" align=middle width=75.52295684999999pt height=30.984656999999984pt/> <img src="/latex/tex/20ae10b1eea36dc5e9febb2da168f429.svg?invert_in_darkmode&sanitize=true" align=middle width=193.59283845pt height=37.80850590000001pt/>.

In the paper we have 
<p align="center"><img src="/latex/tex/b2b4dee3e67cd3cc6f9f2eff2cd9ff17.svg?invert_in_darkmode&sanitize=true" align=middle width=361.7504946pt height=37.663621049999996pt/></p>

* Q3: Most likely I'm mistaken, but I was wondering if (i) in the first fraction the denominator has also the term <img src="/latex/tex/0e1f80bbeed0dd8ac41895e03cf6c1ae.svg?invert_in_darkmode&sanitize=true" align=middle width=42.41839139999999pt height=29.190975000000005pt/> where <img src="/latex/tex/fdba6d2e88a65b633971bf20f13e5034.svg?invert_in_darkmode&sanitize=true" align=middle width=21.00462704999999pt height=24.65753399999998pt/> is the determinant of <img src="/latex/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/>, and (ii) the first <img src="/latex/tex/a781a922ce825b97ac83bb00bf9aee34.svg?invert_in_darkmode&sanitize=true" align=middle width=75.52121444999999pt height=24.65753399999998pt/> is transposed, so that we have
<p align="center"><img src="/latex/tex/bb7f1a5ce6f16824ca87966f1fec6c16.svg?invert_in_darkmode&sanitize=true" align=middle width=415.34634239999997pt height=42.22070985pt/></p>
R3: You are right.


3. **Design/sample**

<img src="/latex/tex/78ec2b7008296ce0561cf83393cb746d.svg?invert_in_darkmode&sanitize=true" align=middle width=14.06623184999999pt height=22.465723500000017pt/> is the random process for the design.

* Q4: D has domain equal to the power set of <img src="/latex/tex/439a14d61ecaa96705b4390d34a96c49.svg?invert_in_darkmode&sanitize=true" align=middle width=21.14196314999999pt height=22.465723500000017pt/>, i.e. <img src="/latex/tex/790d5d68827b3d5165f8914b3650972b.svg?invert_in_darkmode&sanitize=true" align=middle width=50.32573919999999pt height=24.65753399999998pt/>?

R4: we have made the assumption of a fixed size design, so in this case, <img src="/latex/tex/6c0125dedabbb07863101edfd4f2d656.svg?invert_in_darkmode&sanitize=true" align=middle width=37.67357054999999pt height=24.65753399999998pt/> has domain equal to  the power set of <img src="/latex/tex/439a14d61ecaa96705b4390d34a96c49.svg?invert_in_darkmode&sanitize=true" align=middle width=21.14196314999999pt height=22.465723500000017pt/>, i.e. <img src="/latex/tex/790d5d68827b3d5165f8914b3650972b.svg?invert_in_darkmode&sanitize=true" align=middle width=50.32573919999999pt height=24.65753399999998pt/>.
<img src="/latex/tex/78ec2b7008296ce0561cf83393cb746d.svg?invert_in_darkmode&sanitize=true" align=middle width=14.06623184999999pt height=22.465723500000017pt/> is a random variable: <img src="/latex/tex/a3b431d2bd427d589d72aa7d17bc6c75.svg?invert_in_darkmode&sanitize=true" align=middle width=60.64133294999999pt height=22.465723500000017pt/>{the set of all probabilities on <img src="/latex/tex/439a14d61ecaa96705b4390d34a96c49.svg?invert_in_darkmode&sanitize=true" align=middle width=21.14196314999999pt height=22.465723500000017pt/>}, <img src="/latex/tex/c8f0fcddb1c6a7d6921dcc0292039768.svg?invert_in_darkmode&sanitize=true" align=middle width=74.06607779999999pt height=24.65753399999998pt/>


* Q5: the codomain is the set of probability distributions on <img src="/latex/tex/439a14d61ecaa96705b4390d34a96c49.svg?invert_in_darkmode&sanitize=true" align=middle width=21.14196314999999pt height=22.465723500000017pt/>, right? Can we indicate that with <img src="/latex/tex/9de0d1dae382ede944a899902313a9d8.svg?invert_in_darkmode&sanitize=true" align=middle width=28.86213824999999pt height=22.465723500000017pt/>?
R5: the codomain of <img src="/latex/tex/78ec2b7008296ce0561cf83393cb746d.svg?invert_in_darkmode&sanitize=true" align=middle width=14.06623184999999pt height=22.465723500000017pt/> is the set of all probabilities distributions on <img src="/latex/tex/439a14d61ecaa96705b4390d34a96c49.svg?invert_in_darkmode&sanitize=true" align=middle width=21.14196314999999pt height=22.465723500000017pt/>
the codomain of <img src="/latex/tex/6c0125dedabbb07863101edfd4f2d656.svg?invert_in_darkmode&sanitize=true" align=middle width=37.67357054999999pt height=24.65753399999998pt/> is <img src="/latex/tex/acf5ce819219b95070be2dbeb8a671e9.svg?invert_in_darkmode&sanitize=true" align=middle width=32.87674994999999pt height=24.65753399999998pt/>

* Q6: If Q4 and Q5 are right, then we should have <img src="/latex/tex/bbc518087dab66f0f80647ba349dccc6.svg?invert_in_darkmode&sanitize=true" align=middle width=132.52312589999997pt height=24.65753399999998pt/>

R6: <img src="/latex/tex/58cc4acf8355c8bda1a75c81afb528d0.svg?invert_in_darkmode&sanitize=true" align=middle width=140.32528784999997pt height=24.65753399999998pt/>

4. **Exchangeability condition**

* Q7: Why do we use it? (I guess to have property 2.3)
R7: this is a way to ignore the order in which elements are selected. Adaptive sampling is a counter example

* Q8: It is verified in our case? I mean, it's just for SRS or not?
R8: Yes, we limit ourselves to these cases. The selection proportional to size, this is also the case.
