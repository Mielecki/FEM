# funkcja licząca całki kwadraturą gaussa
gaussian_quadrature <- function(f, a, b){
  x1 = -1/sqrt(3)
  x2 = 1/sqrt(3)
  c1 = (b-a)/2
  c2 = (a+b)/2
  return(
    c1*(f(c1*x1+c2)+f(c1*x2+c2))
  )
}

# funkcja zwracająca punkt podziału
x_i <- function(n, i){
  return(2*i/n)
}

# funkcja zwracająca wartość funkcji bazowej w punkcie
e_i <- function(n, i, x){
  if(x_i(n, i-1) <= x && x <= x_i(n, i)){
    return(x*n/2-i+1)
  }
  else if (x_i(n, i) < x && x <= x_i(n, i+1)){
    return(-x*n/2+i+1)
  }
  return(0)
}

# funkcja zwracająca pochodną funkcji bazowej w punkcie
e_i_prim <- function(n, i, x){
  if(x_i(n, i-1) <= x && x <= x_i(n, i)){
    return(n/2)
  }
  else if (x_i(n, i) < x && x < x_i(n, i+1)){
    return(-n/2)
  }
  return(0)
}

# funkcja zwracająca wartość E(x)
e <- function(x){
  if(0 <= x && x <= 1){
    return(3)
  }
  else if(1 < x && x <= 2){
    return(5)
  }
  return(0)
}

# funkcja zwracająca wartość komórki macierzy B
b_i_j <- function(n, i, j){
  a = max(0, x_i(n, i-1), x_i(n, j-1))
  b = min(2, x_i(n, i+1), x_i(n, j+1))
  return(gaussian_quadrature((function(x) e(x)*e_i_prim(n,i,x)*e_i_prim(n, j, x)), a, b) - 3*e_i(n, i, 0)*e_i(n,j,0))
}

# funkcja zwracająca wartość komórki macierzy L
l_i <- function(n, i){
  return(-30*e_i(n,i,0))
}

# funkcja tworząca macierze B i L oraz rozwiązująca je
solution <- function(n){
  B <- matrix(0, nrow=n, ncol=n)
  L <- vector()
  
  for (i in 1:n){
    for (j in 1:n){
      B[i,j] = b_i_j(n, j-1, i-1)
    }
  }

  for (i in 1:n){
    L = c(L, l_i(n, i-1))
  }

  sol = solve(B,L)

  return(
    function(x){
      result =0
      for (i in 1:n){
        result = result + sol[i]*e_i(n,i-1,x)
      }
      return(result)
    }
    )
}

# funkcja rysująca funkcje bazowe
basis_graph <- function(n){
  x = seq(0,2, length.out = 100*n)
  plot(x, mapply(e_i, n, 0, x),
       main = 'funkcje bazowe', type='l', xlab='', ylab='')
  for (i in 1:n){
    lines(x, mapply(e_i, n, i, x))
  }
}

# funkcja rysująca wykres rozwiązania
solution_graph <- function(n){
  x = seq(0,2, length.out = 100*n)
  y = solution(n)
  plot(x, sapply(x, y), main = 'odkształcenie sprężyste', type='l', ylab='',xlab='')
  }

basis_graph(5)
solution_graph(5)