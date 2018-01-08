# Script for plotting nullclines of toggle switch  
# with cooperativity

# Add package for plotting
using Plots
 
k = 100
K = 10 #/s
q = 10
Q = 1
 
f = 0.01
r = 0.01
a = 1
b = 1
 
B=0:0.01:10
 
A1 = k*r/K./(r+f*B.^a)
A2 = (r/f*(q./(Q*B)-1)).^(1/b)

gr() # set plotting back end to gr()

plot(B, [A1 A2], xlabel = "Time", ylabel = "Solution", xlim = (0,10), ylim = (0,10))

savefig("../Results/Nullclines.png")
