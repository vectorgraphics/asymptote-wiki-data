import jacobi;

bool close(pair a, pair b) 
{
	real norm=(b == 0) ? 1 : max(abs(a),abs(b));
	return abs(a-b) <= 1000*realEpsilon*norm;
}

void StartTest(string desc)
{ 
	write("Testing " + desc + "...",flush);
}

void EndTest()
{
	write("PASSED.");
}


StartTest("jacobi package");
assert(close(ljtheta1((2, 2), (0.5, 0.5)), (2.347537035784, -4.509974497964)));
assert(close(ljtheta2((2, 2), (0.5, 0.5)), (2.70710412933166, -2.1537800077716)));
assert(close(ljtheta3((2, 2), (0.5, 0.5)), (2.732256129489196, -2.153780007771601)));
assert(close(ljtheta4((2, 2), (0.5, 0.5)), (2.19990811990507, -6.08077082475884)));
assert(close(kleinj((0, 2)), 66**3));
assert(close(eta((0, 2)), gamma(1/4) / 2**(11/8) / pi**(3/4)));
assert(close(lambda((0.5, 0.5)), 2));
assert(close(ellipticK((0.5, 0.0)), (8 * pi**(3/2) / gamma(-1/4)**2, 0.0)));
assert(close(agm(1, sqrt(2)), 2 * pi**(3/2) * sqrt(2) / gamma(1/4)**2));
assert(close(kleinj(kleinjinv((0.1, 0.2))), (0.1, 0.2)));
assert(close(wp_from_omega((1, 1), (0.25, 0), (0, 2)), (-13.1600232376583, 0)));
EndTest();

