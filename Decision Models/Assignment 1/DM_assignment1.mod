# Set degli indici
set CARGO;
set WAGON;

# Parametri
param profit {CARGO};          # profitto per tonnellata per ciascun tipo di carico
param surface_req {CARGO};     # superficie necessaria per tonnellata per ciascun tipo di carico
param avail {CARGO};           # disponibilità massima (in tonnellate) per ciascun tipo di carico

param weightCap {WAGON};       # capacità massima in peso (tonnellate) per ciascun vagone
param surfaceCap {WAGON};      # capacità massima in superficie per ciascun vagone

# Variabili decisionali:
# x[i,j] rappresenta il numero di tonnellate di carico di tipo i caricate sul vagone j.
var x {CARGO, WAGON} >= 0;

# Funzione obiettivo: massimizzare il profitto totale
maximize TotalProfit:
    sum {i in CARGO, j in WAGON} profit[i] * x[i,j];

# Vincolo di capacità di peso per ogni vagone
subject to WeightCapacity {j in WAGON}:
    sum {i in CARGO} x[i,j] <= weightCap[j];

# Vincolo di capacità di superficie per ogni vagone
subject to SurfaceCapacity {j in WAGON}:
    sum {i in CARGO} surface_req[i] * x[i,j] <= surfaceCap[j];

# Vincolo di disponibilità per ciascun tipo di carico
subject to Availability {i in CARGO}:
    sum {j in WAGON} x[i,j] <= avail[i];

# Comando per risolvere il modello
solve;

# Visualizza le soluzioni ottimali e il profitto totale
display x, TotalProfit;
