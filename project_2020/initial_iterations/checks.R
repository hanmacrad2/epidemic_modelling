
#Check
vecX = vector('numeric',5)
num_days = 5
vecX[1] = 1

for (t in 2:num_days) {
  print(t)
  #Total rate
  print(vecX[1:t-1]) 
  vecX[t] = t
}
