
## New functions or methods

- Implement the `confint` method to compute confidence or credible
  intervals on the model parameters. In the frequentist case, this
  should allow intervals based on profile-likelihood as done in
  **NSGEV**.
  
- Implement `autolayer` for the classes of fitted objects.


## Extend existing functions and methods

- For at least the `poisGPML` case, allow the use of historical data
  following the of **Renext**. The likelihood could include along with
  ordinary `OT` observations: block maxima or largest order statistics
  ("MAX" data), censored-by-level observations ("OTS" data). Mind that
  the terminology used in **Renext** is rather atypical, hence that some
  changes can be required. However the list structure as used in
  ̀Renext::Renouv` to specify blocks seems pretty adequate here.