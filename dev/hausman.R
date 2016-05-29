hausman <- function(fixed,random) {
  #namen der fixed effects im mixed model
  rNames <- names(random@fixef)
  
  #namen der koeffizienten im fixed model
  fNames <- names(coef(fixed))
  
  #koeffizienten die in BEIDEN modellen vorkommen
  timevarNames <- intersect(rNames,fNames)
  
  k <- length(timevarNames)
  
  #varianz matrix aus dem mixed model
  rV <- vcov(random)
  rownames(rV)=rNames
  colnames(rV)=rNames
  
  #differenz der koeffizieten aus dem mixed model - fixed model
  #NUR coefs die in BEIDEN modellen vorkommen (siehe oben)
  bDiff <- (random@fixef)[timevarNames] - coef(fixed)[timevarNames]
  
  #analog für die varianzen
  vDiff <- vcov(fixed)[timevarNames,timevarNames] - rV[timevarNames,timevarNames]
  
  #berechne teststatistik H
  (H <- as.numeric(t(bDiff) %*% solve(vDiff) %*% bDiff))
  
  #ausgabewert teststatistik H und p-wert
  c(H=H,p.value=pchisq((H),k,lower.tail=FALSE))
}

##example
library(lme4)
lmer1 <- lmer(Reaction ~ Days + (1|Subject) , sleepstudy)
lm1 <- lm(Reaction ~ Days + Subject - 1,sleepstudy)

hausman(lm1,lmer1)

require(multcomp)
# gleiches Modell wie lm1
lm1b <- lm(Reaction ~ 0 + as.factor(Subject) +  Days,sleepstudy)

# test mit H_0 fixed effects = random effects
# H_A: fixed effects != random effects
K <- diag(length(coef(lm1b)))[-19,]
rownames(K) <- names(coef(lm1b))[-19]
test <- glht(lm1b, linfct = K, rhs = coef(lmer1)[[1]][,1])
summary(test)
# kein adjustierter p-Wert < 0.05 --> H_0 beibehalten






