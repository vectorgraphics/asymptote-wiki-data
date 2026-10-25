// Copyright 2022 Stéphane Laurent

// Permission is hereby granted, free of charge, to any person obtaining a 
// copy of this software and associated documentation files (the 
// "Software"), to deal in the Software without restriction, including 
// without limitation the rights to use, copy, modify, merge, publish, 
// distribute, sublicense, and/or sell copies of the Software, and to 
// permit persons to whom the Software is furnished to do so, subject to 
// the following conditions:

// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


private real modulo(real a, real p) {
  real i = a > 0 ? floor(a / p) : ceil(a / p);
  return a - i * p;
}

private pair calctheta3(pair z, pair tau) {
  pair out = (1.0, 0.0);
  int n = 0;
  while(true) {
    n = n + 1;
    //real nn = n;
    pair qweight = exp(n * I * pi * (n * tau + 2.0 * z)) +
      exp(n * I * pi * (n * tau - 2.0 * z));
    out += qweight;
    if(length(out) == 0) {
      abort("log(0)");
    } else 
    if(n >= 3 && (out + qweight) == out) {
      break;
    }
  }
  return log(out);
}

private pair argtheta3(pair z, pair tau, int pass_in) {
  int passes = pass_in + 1;
  if(passes > 800) {
    abort("passes > 800 (argtheta3)");
  }
  real z_img = ypart(z);
  real h = ypart(tau) / 2.0;
  pair zuse = (modulo(xpart(z), 1.0), z_img);
  pair out;
  if(z_img < -h) {
    out = argtheta3(-zuse, tau, passes);
  } else if(z_img >= h) {
    pair zmin = zuse - tau;
    out = -2.0 * pi * I * zmin + argtheta3(zmin, tau, passes) -
      I * pi * tau;
  } else {
    out = calctheta3(zuse, tau);
  }
  return out;
}

private pair dologtheta3(pair z, pair tau, int pass_in) {
  int passes = pass_in + 1;
  if(passes > 800) {
    abort("passes > 800 (dologtheta3)");
  }
  pair tau2;
  real rl = xpart(tau);
  if(rl >= 0) {
    tau2 = modulo(rl + 1.0, 2) - 1.0 + I * ypart(tau);
  } else {
    tau2 = modulo(rl - 1.0, 2) + 1.0 + I * ypart(tau);
  }
  rl = xpart(tau2);
  pair out;
  if(rl > 0.6) {
    out = dologtheta3(z + 0.5, tau2 - 1.0, passes);
  } else if(rl < -0.6) {
    out = dologtheta3(z + 0.5, tau2 + 1.0, passes);
  } else if(length(tau2) < 0.98 && ypart(tau2) < 0.98) {
    pair tauprime = (-1.0, 0.0) / tau2;
    out = I * pi * tauprime * z * z +
      dologtheta3(-z * tauprime, tauprime, passes) -
      log(sqrt(-I * tau2));
  } else {
    out = argtheta3(z, tau2, passes);
  }
  return out;
}

private pair M(pair z, pair tau) {
  return I * pi * (z + tau / 4.0);
}

private pair logtheta1(pair z, pair tau) {
  return M(z/pi - 0.5, tau) + dologtheta3(z/pi - 0.5 + 0.5 * tau, tau, 0);
}

/*<asyxml><function type="pair" signature="ljtheta1(pair z, pair tau)"><code></asyxml>*/
pair ljtheta1(pair z, pair tau) {
/*<asyxml></code><documentation>Computes the logarithm of the first Jacobi theta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  if(ypart(tau) <= 0) {
    abort("The complex number `tau` must have a nonnegative imaginary part.");
  }
  return logtheta1(z, tau);
}

private pair logtheta2(pair z, pair tau) {
  return M(z/pi, tau) + dologtheta3(z/pi + 0.5 * tau, tau, 0);
}

/*<asyxml><function type="pair" signature="ljtheta2(pair z, pair tau)"><code></asyxml>*/
pair ljtheta2(pair z, pair tau) {
/*<asyxml></code><documentation>Computes the logarithm of the second Jacobi theta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  if(ypart(tau) <= 0) {
    abort("The complex number `tau` must have a nonnegative imaginary part.");
  }
  return logtheta2(z, tau);
}

private pair logtheta3(pair z, pair tau) {
  return dologtheta3(z/pi, tau, 0);
}

/*<asyxml><function type="pair" signature="ljtheta3(pair z, pair tau)"><code></asyxml>*/
pair ljtheta3(pair z, pair tau) {
/*<asyxml></code><documentation>Computes the logarithm of the third Jacobi theta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  if(ypart(tau) <= 0) {
    abort("The complex number `tau` must have a nonnegative imaginary part");
  }
  return logtheta3(z, tau);
}

private pair logtheta4(pair z, pair tau) {
  return dologtheta3(z/pi + 0.5, tau, 0);
}

/*<asyxml><function type="pair" signature="ljtheta4(pair z, pair tau)"><code></asyxml>*/
pair ljtheta4(pair z, pair tau) {
/*<asyxml></code><documentation>Computes the logarithm of the fourth Jacobi theta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  if(ypart(tau) <= 0) {
    abort("The complex number `tau` must have a nonnegative imaginary part.");
  }
  return logtheta4(z, tau);
}

/*<asyxml><function type="pair" signature="jtheta1(pair z, pair tau)"><code></asyxml>*/
pair jtheta1(pair z, pair tau) {
/*<asyxml></code><documentation>Computes the first Jacobi theta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  return exp(ljtheta1(z, tau));
}

/*<asyxml><function type="pair" signature="jtheta2(pair z, pair tau)"><code></asyxml>*/
pair jtheta2(pair z, pair tau) {
/*<asyxml></code><documentation>Computes the second Jacobi theta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  return exp(ljtheta2(z, tau));
}

/*<asyxml><function type="pair" signature="jtheta3(pair z, pair tau)"><code></asyxml>*/
pair jtheta3(pair z, pair tau) {
/*<asyxml></code><documentation>Computes the third Jacobi theta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  return exp(ljtheta3(z, tau));
}

/*<asyxml><function type="pair" signature="jtheta4(pair z, pair tau)"><code></asyxml>*/
pair jtheta4(pair z, pair tau) {
/*<asyxml></code><documentation>Computes the fourth Jacobi theta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  return exp(ljtheta4(z, tau));
}

private pair _eta_(pair tau) {
  return exp(I * pi * tau/12.0 + dologtheta3((tau + 1.0)/2.0, 3.0 * tau, 0)); 
}

/*<asyxml><function type="pair" signature="eta(pair tau)"><code></asyxml>*/
pair eta(pair tau) {
/*<asyxml></code><documentation>Computes the Dedekind eta function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  if(ypart(tau) <= 0) {
    abort("The complex number `tau` must have a nonnegative imaginary part.");
  }
  return _eta_(tau);
}

private pair theta1(pair z, pair tau) {
  return exp(logtheta1(z, tau));
}

private pair theta2(pair z, pair tau) {
  return exp(logtheta2(z, tau));
}

private pair theta3(pair z, pair tau) {
  return exp(logtheta3(z, tau));
}

private pair theta4(pair z, pair tau) {
  return exp(logtheta4(z, tau));
}

/*<asyxml><function type="pair" signature="lambda(pair tau)"><code></asyxml>*/
pair lambda(pair tau) {
/*<asyxml></code><documentation>Computes the modular lambda function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  if(ypart(tau) <= 0) {
    abort("The complex number `tau` must have a nonnegative imaginary part.");
  }
  return (theta2(0, tau) / theta3(0, tau))**4; 
}

/*<asyxml><function type="pair" signature="kleinj(pair tau)"><code></asyxml>*/
pair kleinj(pair tau) {
/*<asyxml></code><documentation>Computes the Klein j-invariant function.
  The imaginary part of 'tau' must be nonnegative.</documentation></function></asyxml>*/
  if(ypart(tau) <= 0) {
    abort("The complex number `tau` must have a nonnegative imaginary part.");
  }
  pair lbd = lambda(tau);
  pair x = lbd * (1.0 - lbd);
  return 256.0 * (1.0 - x)**3 / x**2; 
}

/*<asyxml><function type="pair" signature="CarlsonRF(pair x, pair y, pair z)"><code></asyxml>*/
pair CarlsonRF(pair x, pair y, pair z) {
/*<asyxml></code><documentation>Computes the Carlson integral 'RF'.
  At most one of 'x', 'y', 'z' can be zero.</documentation></function></asyxml>*/
  int xzero = x == (0, 0) ? 1 : 0; 
  int yzero = y == (0, 0) ? 1 : 0; 
  int zzero = z == (0, 0) ? 1 : 0; 
  int nzeros = xzero + yzero + zzero;
  if(nzeros >= 2) {
    abort("At most one of `x`, `y`, `z` can be 0.");
  }
  real dx = realMax;
  real dy = realMax;
  real dz = realMax;
  real epsilon = 2.0 * realEpsilon;
  pair A;
  while(dx > epsilon || dy > epsilon || dz > epsilon) {
    pair lambda = sqrt(x)*sqrt(y) + sqrt(y)*sqrt(z) + sqrt(z)*sqrt(x);
    x = (x + lambda) / 4.0;
    y = (y + lambda) / 4.0;
    z = (z + lambda) / 4.0;
    A = (x + y + z) / 3.0;
    dx = length(1.0 - x/A);
    dy = length(1.0 - y/A);
    dz = length(1.0 - z/A);
  } 
  real E2 = dx*dy + dy*dz + dz*dx;
  real E3 = dy*dx*dz;
  return (1 - E2/10 + E3/14 + E2*E2/24 - 3*E2*E3/44 - 5*E2*E2*E2/208 +
    3*E3*E3/104 + E2*E2*E3/16) / sqrt(A);
}

/*<asyxml><function type="pair" signature="CarlsonRD(pair x, pair y, pair z)"><code></asyxml>*/
pair CarlsonRD(pair x, pair y, pair z) {
/*<asyxml></code><documentation>Computes the Carlson integral 'RD'.
  At most one of 'x', 'y', 'z' can be zero.</documentation></function></asyxml>*/
  int xzero = x == (0, 0) ? 1 : 0; 
  int yzero = y == (0, 0) ? 1 : 0; 
  int zzero = z == (0, 0) ? 1 : 0; 
  int nzeros = xzero + yzero + zzero;
  if(nzeros >= 2) {
    abort("At most one of `x`, `y`, `z` can be 0.");
  }
  real dx = realMax;
  real dy = realMax;
  real dz = realMax;
  real epsilon = 2.0 * realEpsilon;
  pair s = (0.0, 0.0);
  pair fac = (1.0, 0.0);
  pair A;
  while(dx > epsilon || dy > epsilon || dz > epsilon) {
    pair lambda = sqrt(x)*sqrt(y) + sqrt(y)*sqrt(z) + sqrt(z)*sqrt(x);
    s = s + fac/(sqrt(z) * (z + lambda));
    fac = fac / 4.0;
    x = (x + lambda) / 4.0;
    y = (y + lambda) / 4.0;
    z = (z + lambda) / 4.0;
    A = (x + y + 3.0 * z) / 5.0;
    dx = length(1.0 - x/A);
    dy = length(1.0 - y/A);
    dz = length(1.0 - z/A);
  } 
  real E2 = dx * dy + dy * dz + 3 * dz * dz + 2 * dz * dx +  
              dx * dz + 2 * dy * dz;
  real E3 = dz * dz * dz + dx * dz * dz + 3 * dx * dy * dz + 
              2 * dy * dz * dz + dy * dz * dz + 2 * dx * dz * dz;
  real E4 = dy * dz * dz * dz + dx * dz * dz * dz + dx * dy * dz * dz + 
              2 * dx * dy * dz * dz;
  real E5 = dx * dy * dz * dz * dz;
  return 3 * s + fac * (1 - 3 * E2/14 + E3/6 + 9 * E2 * E2/88 - 3 * E4/22 - 
          9 * E2 * E3/52 + 3 * E5/26 - E2 * E2 * E2/16 + 
          3 * E3 * E3/40 + 3 * E2 * E4/20 + 45 * E2 * E2 * E3/272 - 
          9 * (E3 * E4 + E2 * E5)/68) / A / sqrt(A);
}

private pair atanhz(pair z) {
  return 0.5 * log((1+z) / (1-z));
} 

private bool isinf(pair z) {
  return abs(xpart(z)) == inf || abs(ypart(z)) == inf;
} 

private real rfloor(real x) {
	return x - (x % 1);
}

private real rceil(real x) {
	return rfloor(x) + 1;
}

/*<asyxml><function type="pair" signature="ellipticF(pair phi, pair m)"><code></asyxml>*/
pair ellipticF(pair phi, pair m) {
/*<asyxml></code><documentation>Computes the incomplete elliptic integral of first kind; 
  'phi' is the amplitude and 'm' is the squared modulus.</documentation></function></asyxml>*/
  if(phi == 0 || m == inf || m == -inf){
    return (0, 0);
  }
  if(xpart(phi) == 0 && ypart(phi) == inf && ypart(m) == 0 && xpart(m) > 0 &&
      xpart(m) < 1){
    return sgn(ypart(phi)) *
      (ellipticF(pi/2, m) - ellipticF(pi/2, 1/m) / sqrt(m));
  }
  if(abs(xpart(phi)) == pi/2 && m == 1){
    return (nan, nan);
  }
  if(xpart(phi) >= -pi/2 && xpart(phi) <= pi/2){
    if(m == 1 && abs(xpart(phi)) < pi/2){
      return atanhz(sin(phi)); 
    }
    if(m == 0){
      return phi;
    }
    pair sine = sin(phi);
    if(isinf(sine)) {
      abort("`sin(phi)` is not finite.");
    }
    pair sine2 = sine*sine;
    pair cosine2 = 1 - sine2;
    pair oneminusmsine2 = 1 - m * sine2;
    return sine * CarlsonRF(cosine2, oneminusmsine2, (1.0, 0.0));
  }
  real k;
  if(xpart(phi) > pi/2) {
    k = rceil((xpart(phi)-pi/2) / pi);
    phi = phi - k * pi;
    // while(xpart(phi) > pi/2){
    //   phi = phi - pi;
    //   k = k + 1;
    // }
  } else {
    k = -rfloor((pi/2-xpart(phi)) / pi);
    phi = phi - k * pi;
    // while(xpart(phi) < -pi/2){
    //   phi = phi + pi;
    //   k = k - 1;
    // }
  }
  return 2*k*ellipticF(pi/2, m) + ellipticF(phi, m);
}

/*<asyxml><function type="pair" signature="ellipticK(pair m)"><code></asyxml>*/
pair ellipticK(pair m) {
/*<asyxml></code><documentation>Computes the complete elliptic integral of first kind; 
  'm' is the squared modulus.</documentation></function></asyxml>*/
  return ellipticF((pi/2, 0.0), m);
}

/*<asyxml><function type="pair" signature="ellipticE(pair phi, pair m)"><code></asyxml>*/
pair ellipticE(pair phi, pair m) {
/*<asyxml></code><documentation>Computes the incomplete elliptic integral of second kind; 
  'phi' is the amplitude and 'm' is the squared modulus.</documentation></function></asyxml>*/
  if(phi == 0){
    return (0, 0);
  }
  if(xpart(m) == inf && ypart(m) == 0){
    return (nan, nan);
  }
  if(xpart(phi) >= -pi/2 && xpart(phi) <= pi/2){
    if(m == (0, 0)){
      return phi;
    }
    if(m == (1, 0)){
      return sin(phi);
    }
    pair sine = sin(phi);
    if(isinf(sine)){
      abort("`sin(phi)` is not finite.");
    }
    pair sine2 = sine*sine;
    pair cosine2 = 1 - sine2;
    pair oneminusmsine2 = 1 - m * sine2;
    return sine * (CarlsonRF(cosine2, oneminusmsine2, (1, 0)) -
            m * sine2 * CarlsonRD(cosine2, oneminusmsine2, (1, 0)) / 3);
  }
  real k;
  if(xpart(phi) > pi/2) {
    k = rceil((xpart(phi)-pi/2) / pi);
    phi = phi - k * pi;
    // while(xpart(phi) > pi/2){
    //   phi = phi - pi;
    //   k = k + 1;
    // }
  } else {
    k = -rfloor((pi/2-xpart(phi)) / pi);
    phi = phi - k * pi;
    // while(xpart(phi) < -pi/2){
    //   phi = phi + pi;
    //   k = k - 1;
    // }
  }
  return 2*k*ellipticE(pi/2, m) + ellipticE(phi, m);
}

/*<asyxml><function type="pair" signature="ellipticE(pair m)"><code></asyxml>*/
pair ellipticE(pair m) {
/*<asyxml></code><documentation>Computes the complete elliptic integral of second kind; 
  'm' is the squared modulus.</documentation></function></asyxml>*/
  return ellipticE((pi/2, 0.0), m);
}

private pair Eisen2(pair tau) {
  pair j3sq = theta3(0, tau)**2;
  return 6/pi * ellipticE(lambda(tau)) * j3sq - j3sq**2 - theta4(0, tau)**4;
} 

/*<asyxml><function type="pair" signature="EisensteinE2(pair q)"><code></asyxml>*/
pair EisensteinE2(pair q) {
/*<asyxml></code><documentation>Computes the Eisenstein E-series of weight 2.
  The modulus of 'q' must be strictly smaller than 1, and 'q' must not be a negative real number.</documentation></function></asyxml>*/
  if(length(q) >= 1.0) {
    abort("The modulus of `q` must be strictly smaller than 1.");
  }
  if(ypart(q) == 0.0 && xpart(q) <= 0) {
    abort("If `q` is real, it must be strictly positive.");
  }
  pair tau = -I * log(q) / pi / 2.0;
  return Eisen2(tau); 
}

/*<asyxml><function type="pair" signature="EisensteinE4(pair q)"><code></asyxml>*/
pair EisensteinE4(pair q) {
/*<asyxml></code><documentation>Computes the Eisenstein E-series of weight 4.
  The modulus of 'q' must be strictly smaller than 1, and 'q' must not be a negative real number.</documentation></function></asyxml>*/
  if(length(q) >= 1.0) {
    abort("The modulus of `q` must be strictly smaller than 1.");
  }
  if(ypart(q) == 0.0 && xpart(q) <= 0) {
    abort("If `q` is real, it must be strictly positive.");
  }
  pair tau = -I * log(q) / pi / 2.0;
  return (theta2(0, tau)**8 + theta3(0, tau)**8 + theta4(0, tau)**8) / 2.0; 
}

/*<asyxml><function type="pair" signature="EisensteinE6(pair q)"><code></asyxml>*/
pair EisensteinE6(pair q) {
/*<asyxml></code><documentation>Computes the Eisenstein E-series of weight 6.
  The modulus of 'q' must be strictly smaller than 1, and 'q' must not be a negative real number.</documentation></function></asyxml>*/
  if(length(q) >= 1.0) {
    abort("The modulus of `q` must be strictly smaller than 1.");
  }
  if(ypart(q) == 0.0 && xpart(q) <= 0) {
    abort("If `q` is real, it must be strictly positive.");
  }
  pair tau = -I * log(q) / pi / 2.0;
  return (theta3(0, tau)**12 + theta4(0, tau)**12 - 3.0 * theta2(0, tau)**8 * 
    (theta3(0, tau)**4 + theta4(0, tau)**4)) / 2.0; 
}

/*<asyxml><function type="pair" signature="agm(pair x, pair y)"><code></asyxml>*/
pair agm(pair x, pair y) {
/*<asyxml></code><documentation>Computes the arithmtic-geomeric mean.</documentation></function></asyxml>*/
  if(x + y == 0 || x == 0 || y == 0){
    return((0, 0));
  }
  return pi/4 * (x + y) / ellipticK(((x-y)/(x+y))**2);
}

/*<asyxml><function type="pair" signature="kleinjinv(pair j)"><code></asyxml>*/
pair kleinjinv(pair j) {
/*<asyxml></code><documentation>Computes the inverse of the Klein j-invariant function.</documentation></function></asyxml>*/
  pair x;
  if(isinf(j)) {
    x = (0, 0);
  } else {
    pair t = (-j*j*j + 2304 * j*j + 12288 *
          sqrt(3 * (1728 * j*j - j*j*j)) - 884736 * j)**(1/3);
    x = 1/768 * t - (1536 * j - j*j) / (768 * t) + (1 - j/768);
  }
  pair lbd = -(-1 - sqrt(1 - 4*x)) / 2;
  return I * agm(1, sqrt(1-lbd)) / agm(1, sqrt(lbd));
}

private bool areclose(pair z1, pair z2) {
  real lz2 = length(z2);
  real norm = (lz2 < realEpsilon) ? 1.0 : max(length(z1), lz2);
  return length(z1 - z2) < 2.0 * realEpsilon * norm;
}

private pair ipower(pair z, int p) {
  if(p == 0) {
    return (1.0, 0.0);
  }
  if(p == 1) {
    return z;
  }
  pair tmp = ipower(z, p # 2);
  if(p % 2 == 0) {
    return tmp * tmp;
  } else {
    return z * tmp * tmp;
  }
}

private pair jtheta1dash(pair z, pair tau) {
  pair q = exp(I * pi * tau);
  pair out = (0.0, 0.0);
  pair alt = (-1.0, 0.0);
  for(int n = 0; n < 1000; ++n) {
    alt = -alt;
    int k = 2 * n + 1;
    pair outnew = out + alt * ipower(q, n * (n + 1)) * k * cos(k * z);
    if(areclose(out, outnew)) {
      return 2.0 * q**0.25 * out;
    }
    out = outnew;
  }
  abort("Reached 1000 iterations (theta1dash).");
  return (0.0, 0.0);
}

private pair jtheta1dash0(pair tau){
  return exp(logtheta2((0.0, 0.0), tau) + logtheta3((0.0, 0.0), tau) + 
    logtheta4((0.0, 0.0), tau));
}

private pair jtheta1dashdashdash0(pair tau){
  return -2.0 * _eta_(tau)**3 * Eisen2(tau);
}

private pair dljtheta1(pair z, pair tau) {
  if(z == (0, 0)){
    return(jtheta1dash0(tau) / theta1((0,0), tau));
  }
  return jtheta1dash(z, tau) / theta1(z, tau);
}

private pair wp_from_tau(pair z, pair tau) { // wp(z, omega1 = 1/2, omega2 = tau/2)
  pair j2 = theta2(0, tau);
  pair j3 = theta3(0, tau);
  pair j1 = theta1(z*pi, tau);
  pair j4 = theta4(z*pi, tau);
  return (pi * j2 * j3 * j4 / j1)**2 - pi**2 * (j2**4 + j3**4) / 3.0;
}

pair wp_from_omega(pair z, pair omega1, pair omega2) {
  return wp_from_tau(z/omega1/2, omega2/omega1) / omega1 / omega1 / 4;
}