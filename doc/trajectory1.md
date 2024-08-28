Trajectory modeling steps 2-4: representation, scoring, and search process {#trajectory1}
====================================

Now, we have snapshot models for various intermediate states along our process of interest. In this step, we will combine these snapshot models to produce trajectories.

# Representing the model

# Scoring the model

# Searching for good scoring models

 \f[
 W(\mathcal{X}) \propto   \displaystyle&\prod^{T}_{t=0} \mathcal{P}( X_{N,t}, N_{t} | D_{t}) \cdot\\ \displaystyle&\prod^{T-1}_{t=0} W(X_{N,t+1},N_{t+1} | X_{N,t},N_{t}, D_{t,t+1}),
 \f]