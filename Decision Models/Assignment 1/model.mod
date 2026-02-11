# set of the indexes
set CARGO;
set WAGON;

# parameters
param profit {CARGO};        
param surface_req {CARGO};    
param avail {CARGO};          

param weightCap {WAGON};     
param surfaceCap {WAGON};      

# decision variables:
var x {CARGO, WAGON} >= 0;

# objective function (maximize total profit)
maximize TotalProfit:
    sum {i in CARGO, j in WAGON} profit[i] * x[i,j];


# constraints 

# weight capacity constraint
subject to WeightCapacity {j in WAGON}:
    sum {i in CARGO} x[i,j] <= weightCap[j];

# surface capacity constraint
subject to SurfaceCapacity {j in WAGON}:
    sum {i in CARGO} surface_req[i] * x[i,j] <= surfaceCap[j];

# availability constraint
subject to Availability {i in CARGO}:
    sum {j in WAGON} x[i,j] <= avail[i];

# solve the model and visualize the optimal solution
solve;

display x, TotalProfit;
