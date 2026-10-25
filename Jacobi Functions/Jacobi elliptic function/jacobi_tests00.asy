
import plain_pens;
import palette;
import jacobi;

size(1024, 1024);

real modulo(real a, real p) {
  real i = a > 0 ? floor(a / p) : ceil(a / p);
  return a - i * p;
}

pair calctheta3(pair z, pair tau) {
	pair out = (1.0, 0.0);
	int n = 0;
	while(true) {
		n = n + 1;
		real nn = n;
		pair qweight = exp(nn * I * pi * (nn * tau + 2.0 * z)) +
			exp(nn * I * pi * (nn * tau - 2.0 * z));
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

pair argtheta3(pair z, pair tau, int pass_in) {
	int passes = pass_in + 1;
	if(passes > 800) {
		abort("passes > 800 (argtheta3)");
	}
	real z_img = z.y;
	real h = tau.y / 2.0;
	pair zuse = (modulo(z.x, 1.0), z_img);
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

pair dologtheta3(pair z, pair tau, int pass_in) {
	int passes = pass_in + 1;
	if(passes > 800) {
		abort("passes > 800 (dologtheta3)");
	}
	pair tau2;
	real rl = tau.x;
	if(rl >= 0) {
		tau2 = modulo(rl + 1.0, 2) - 1.0 + I * tau.y;
	} else {
		tau2 = modulo(rl - 1.0, 2) + 1.0 + I * tau.y;
	}
	rl = tau2.x;
	pair out;
	if(rl > 0.6) {
		out = dologtheta3(z + 0.5, tau2 - 1.0, passes);
	} else if(rl < -0.6) {
		out = dologtheta3(z + 0.5, tau2 + 1.0, passes);
	} else if(length(tau2) < 0.98 && tau2.y < 0.98) {
		pair tauprime = (-1.0, 0.0) / tau2;
		out = I * pi * tauprime * z * z +
			dologtheta3(-z * tauprime, tauprime, passes) -
			log(sqrt(-I * tau2));
	} else {
		out = argtheta3(z, tau2, passes);
	}
	return out;
}

pair M(pair z, pair tau) {
	return I * pi * (z + tau / 4.0);
}

pair ljtheta1(pair z, pair tau) {
	return M(z/pi - 0.5, tau) + dologtheta3(z/pi -0.5 + 0.5 * tau, tau, 0);
}

pair ljtheta2(pair z, pair tau) {
	return M(z/pi, tau) + dologtheta3(z/pi + 0.5 * tau, tau, 0);
}

pair ljtheta3(pair z, pair tau) {
	return dologtheta3(z/pi, tau, 0);
}

pair ljtheta4(pair z, pair tau) {
	return dologtheta3(z/pi + 0.5, tau, 0);
}

pair jtheta1(pair z, pair tau) {
	return exp(ljtheta1(z, tau));
}

pair jtheta2(pair z, pair tau) {
	return exp(ljtheta2(z, tau));
}

pair jtheta3(pair z, pair tau) {
	return exp(ljtheta3(z, tau));
}

pair jtheta4(pair z, pair tau) {
	return exp(ljtheta4(z, tau));
}

pair eta(pair tau) {
	return exp(I * pi * tau/12.0 + dologtheta3((tau + 1.0)/2.0, 3.0 * tau, 0));
}

// checks:
write(ljtheta1((2, 2), (0.5, 0.5))); // 2.347537-4.509974i
write(ljtheta2((2, 2), (0.5, 0.5))); // 2.707104-2.15378i
write(ljtheta3((2, 2), (0.5, 0.5))); // 2.732256-2.15378i
write(ljtheta4((2, 2), (0.5, 0.5))); // 2.199908-6.080771i
write(eta((0.5, 0.5))); // 0.9057633+0.119246i

pen colormap1(pair z) {
	real x = z.x;
	real y = z.y;
	real modulus = length(z);
	real a = angle(z, warn=false);
	real r = modulo(modulus, 1.0);
	real g = abs(modulo(a, 0.5)) * 2.0;
	real b = abs(modulo(x * y, 1.0));
	return rgb(
		(1.0 - cos(r - 0.5)) * 8.0,
		(1.0 - cos(g - 0.5)) * 8.0,
		(1.0 - cos(b - 0.5)) * 8.0
	);
}

pen colormap2(pair z) {
	real modulus = length(z);
	real h = degrees(z, warn=false);
	real s = sqrt((1+sin(2*pi*log1p(modulus)))/2);
	real v = (1+cos(2*pi*log1p(modulus)))/2;
	return hsv(h,s,v);
}

int n = 128;
pen[][] v=new pen[n][n];
pair a = (-1.0, -1.0);
pair b = (1.0, 1.0);
real dx = (b-a).x / n;
real dy = (b-a).y / n;

for(int i = 0; i < n; ++i){
	for(int j = 0; j < n; ++j){
		pair q = a + (i * dx, j * dy);
		if(length(q) > 0.99999 || (q.y == 0.0 && q.x <= 0.0)) {
			v[i][j] = rgb("15191eff");
		} else {
			pair tau = -I * log(q) / pi;
			v[i][j] = colormap1(eta(tau));
		}
	}
}

image(v, (0, 0), (1, 1));
