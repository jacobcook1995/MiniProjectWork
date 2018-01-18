# Script for plotting nullclines of toggle switch
# with cooperativity

# Add package for plotting
using Plots

const k = 100
const K = 10
const q = 10
const Q = 1
const r = 100000
const f = 10000
a = 2
b = 2

B = 0:0.01:10

A1 = k*r/K./(r+f*B.^a)
A2 = (r/f*(q./(Q*B)-1)).^(1/b)

gr() # set plotting back end to gr()

plot(B, [A1 A2], xlabel = "Time", ylabel = "Solution", xlim = (0,10), ylim = (0,10))

savefig("../Results/Nullclines.png")
