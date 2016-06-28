function res=bisection(fun,a,b,tol,params)

  if fun(a,params)*fun(b,params)  > 0
    error('function has same signs at both endpoints')
  endif
  
  while b-a > tol
    x = (a+b)/2;
    y = fun(x,params);

    if y == 0 
      res = x;
      return
    endif

    if fun(a,params)*y < 0
      b = x;
    else
      a = x;
    endif

  endwhile
  
  res= (b+a)/2;

endfunction
