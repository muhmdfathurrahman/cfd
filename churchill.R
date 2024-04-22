churchill <- function(Re, epsilon, D) { 
  # A function to calculate Darcy friction factor f as a function of Reynolds
  # Number Re.
  
  A = (2.457 * log((7/Re)^.9 + .27 * ( epsilon/D)))^16
  B = (37530 / Re)^16
  return(8 * ((8 / Re)^12 + (A + B)^(-1.5))^(1 / 12))
}


calcRe <- function(Speed, HydraulicDiameter, DynamicViscosity) {
  return(Speed * HydraulicDiameter / DynamicViscosity)
}


U_func1 <- function(H, f, L, D, KL, g = 9.81) {
  return(sqrt(2 * g * H / (f * (L/D) + KL)))
}


churchill_converge1 <- function(firstguess, nu, epsilon, H, L, D, KL) {
  # Converge f by calculating U using the U_func1 function
  oldf = firstguess
  i = 1
  repeat {
    U = U_func1(H, oldf, L, D, KL)
    Re = calcRe(U, D, nu)
    newf = churchill(Re, epsilon, D)
    print(paste("i    :", i))
    print(paste("oldf :", oldf))
    print(paste("U    :", U))
    print(paste("Re   :", Re))
    print(paste("newf :", newf))
    print("=======================")
    if(oldf == newf) {break} else {oldf = newf; i = i + 1}
  }
  return(newf)
}