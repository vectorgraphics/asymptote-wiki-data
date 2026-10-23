
import plain_pens;
import palette;
import jacobi;
size(2048, 2048);

real rfloor(real x) {
	return x - (x % 1);
}

real rceil(real x) {
	return rfloor(x) + 1;
}

real modulo(real a, real p) {
	real i = a > 0 ? rfloor(a / p) : rceil(a / p);
	return a - i * p;
}

pen colormap1(pair z) {
	real x = z.x;
	real y = z.y;
	real modulus = length(z);
	real a = angle(z, warn=false); 
	real r = modulus % 1.0;
	real g = abs(modulo(a, 0.5)) * 2.0;
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
	real s = sqrt((1+sin(2*pi*log1p(modulus)))/2);
	real v = (1+cos(2*pi*log1p(modulus)))/2;
	return hsv(h,s,v);
}

int n = 1024;
pen[][] v=new pen[n][n];
pair a = (-1.0, -1.0);
pair b = (1.0, 1.0);
real dx = (b-a).x / n;
real dy = (b-a).y / n;

for(int i = 0; i < n; ++i){
	for(int j = 0; j < n; ++j){
		pair z = a + (i * dx, j * dy);
		v[i][j] = colormap1(wp_from_omega(z, (0.5, 0), (0.0, 0.5)));
	}
}

image(v, (0, 0), (1, 1));
