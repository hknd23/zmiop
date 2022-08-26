# `zmiop`
*Benjamin E. Bagozzi, Minnie M. Joo, Nguyen K Huynh, and Bumba Mukherjee*

An excessive (“inflated”) share of observations—stemming from two distinct d.g.p’s—fall into a single choice category in many ordered and unordered polytomous outcome variables. Standard Ordered Probit models cannot account for such category inflation which leads to biased inferences. The package offers tools to:


Fit the Zero-Inflated Ordered Probit (ZIOP) model without and with correlated errors (ZIOPC model) to evaluate zero-inflated ordered choice outcomes that result from a dual data generating process (d.g.p.).

Fit the Middle-Inflated Ordered Probit (MIOP) model without and with correlated errors (MIOPC) to account for the inflated middle-category in ordered choice measures related to a dual d.g.p.

### Installation

```
devtools::install_github("hknd23/zmiop")
library(zmiop)
```

### ZiOP Model
```
bp   <- BP_AER_2009
#bp   <- BP_AER_2009[sample(1:nrow(BP_AER_2009), 500, replace = FALSE), ]
bp   <- na.omit(bp)
vars <- c('logGDPpc', 'parliament', 'disaster', 'major_oil', 'major_primary')
Y    <- bp$rep_civwar_DV
X    <- as.matrix(bp[vars])
Z    <- as.matrix(bp[vars])
```
```
## ZIOP
model1_ziop <- iop(Y ~ X | Z, data = bp[vars], type = c("ziop"))
summary(model1_ziop)
model_ziop$coefficients
```

### MiOP Model

```
EUK  <- EUKnowledge
#EUK  <- EUKnowledge[sample(1:nrow(EUKnowledge), 1000, replace = FALSE), ]
EUK  <- na.omit(EUK)
vars <- c('polit_trust', 'Xenophobia', 'discuss_politics', 'univers_ed', 'Professional',
          'Executive', 'Manual', 'Farmer', 'Unemployed', 'rural', 'female', 'age',
          'student', 'income', 'Educ_high', 'Educ_high_mid', 'Educ_low_mid')
Y    <- EUK$EU_support
X    <- as.matrix(EUK[vars])
Z    <- as.matrix(EUK[vars])
```


```
model_miop <- iop(Y ~ X | Z, data = EUK[vars], type = c("miop"))
summary(model_miop)
model_miop$coefficients
```
