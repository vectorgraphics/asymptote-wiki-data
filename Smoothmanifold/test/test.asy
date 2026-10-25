size(7cm);


// pair cp = (pos + 2 * radius * dir(angle));
// pair cp = (2.18299236741329,0.5849310419317);

// pair cp = (2.182992367,0.584931041);
pair cp = (2.18299236,0.58493104);

// pair cp = (2.1829923,0.5849310);
// pair cp = (2.182992,0.584931);
// pair cp = (2.1, 0.5);

path gs = (0,0){(1,0)} .. cp .. {(0,-1)}cycle;

draw(gs);
dot(cp, blue+2);
