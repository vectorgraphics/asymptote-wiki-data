import plain_pens;
import palette;
import "Jacobi Functions/Jacobi theta functions and Weierstrass p/jacobi.asy" as jacobi;
size(1024, 1024);

pen colormap1(pair z) {
	real x = z.x;
	real y = z.y;
	real modulus = length(z);
	real a = angle(z, warn=false);
	real r = modulus % 1.0;
	real g = (a % 0.5) * 2.0;
	real b = (x * y) % 1.0;
	return rgb(
		(1.0 - cos(r - 0.5)) * 8.0,
		(1.0 - cos(g - 0.5)) * 8.0,
		(1.0 - cos(b - 0.5)) * 8.0
	);
}

pen colormap2(pair z) {
	real modulus = length(z);
	real h = degrees(z, warn=false);
	real x = 2*pi*log1p(modulus);
	real s = sqrt((1+sin(x))/2);
	real v = (1+cos(x))/2;
	return hsv(h,s,v);
}

int n = 128;
pen[][] v=new pen[n][n];
pair a = (-1.0, -1.0);
pair b = (1.0, 1.0);
real dx = (b - a).x / n;
real dy = (b - a).y / n;

for(int i = 0; i < n; ++i){
	for(int j = 0; j < n; ++j){
		pair q = a + (i * dx, j * dy);
		if(length(q) > 0.99999 || (q.y == 0.0 && q.x <= 0.0)) {
			v[i][j] = rgb("15191eff");
		} else {
			pair tau = -I * log(q) / pi;
			v[i][j] = colormap1(lambda(tau));
		}
	}
}

image(v, (0, 0), (1, 1));