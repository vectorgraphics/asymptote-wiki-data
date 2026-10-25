size(7cm);


// pair cp = (pos + 2 * radius * dir(angle));
// pair cp = (2.18299236741329,0.5849310419317);
pair cp = (2.182992367,0.584931041);
// pair cp = (2.182992,0.584931);
// pair cp = (2.1, 0.5);

dot(cp, blue+2);

write(cp);

// path gs = pos{dir1} .. (pos + 2 * radius * (1+.1*unitrand())*dir(angle)) .. {dir2}pos;
path gs = (0,0){(1,0)} .. cp .. {(0,-1)}cycle;

draw(gs, red);
