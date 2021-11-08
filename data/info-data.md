# info-data

## Scenarios
1. 0 prognostiche e 10 predittive (come confronto con Ma et al. 2016, Biometrics dove il modello prevede solo markers predittivi)

2. 2 prognostiche (nella scala originale) e 10 predittive. The file is called `scenario2.rda`.
  
  - 2 mod 2 prognostiche (nella scala originale) e 20 predittive. 
    it is stored in `modscenario2.rda`

3. 2 prognostiche (with power trasformation) e 10 predittive. Da quello che ho capito in Ma et al. 2019 BiomJ questo scenario è confrontato con il b) per verificare che un effetto più forte delle covariate prognostiche impatti il Relative Gain in Treatment utility with respect to the other treatment (%MTUg) senza necessariamente alterare la direzione del benefit (MOT). Quindi credo che rispetto al b) vorremmo vedere MOT simile e %MTUg più grande.

4. 2 prognostiche (nella scala originale) e 25 predittive. Questo è un tentativo per aumentare il numero di marker predittivi.

5. 2 prognostiche (nella scala originale) e 25 predittive, di cui 10 effettivamente usate per generare la risposta e le restanti di rumore. Questo è un tentativo per aumentare il numero di marker predittivi.