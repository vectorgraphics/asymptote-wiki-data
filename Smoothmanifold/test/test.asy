size(7cm);


// srand(0);

// pair randomdir (pair dir, real angle)
// { return dir(degrees(dir) + (unitrand()-.5)*angle); }

// path randompath (pair[] controlpoints, real angle)
// // Like `connect`, but the path has some randomness to it.
// {
//     if (controlpoints.length < 2) return nullpath;
//
//     pair outdir = randomdir(controlpoints[1]-controlpoints[0], angle);
//     guide res = controlpoints[0];
//
//     for (int i = 1; i < controlpoints.length-1; ++i)
//     {
//         pair indir = randomdir(controlpoints[(i+1)]-controlpoints[i-1], angle);
//         res = res{outdir} .. {indir}controlpoints[i];
//         outdir = indir;
//     }
//
//     pair indir = randomdir(controlpoints[(controlpoints.length-1)]-controlpoints[controlpoints.length-2], angle);
//
//     return (path)(res{outdir}..{indir}controlpoints[controlpoints.length-1]);
// }

// pair pos = (1,1);

// real angle = 15;
// real radius = 1.13;
// real range = 30;

// pair dir1 = (1,-.5);
// pair dir2 = (-.3,-1);

// draw(pos -- (pos + dir1), purple);
// draw(pos -- (pos + dir2), green);

// pair cp = (pos + 2 * radius * dir(angle));
pair cp = (3.18299236741329,1.5849310419317);
// pair cp = (3.1, 1.5);

dot(cp, blue+2);

write(cp);

// path gs = pos{dir1} .. (pos + 2 * radius * (1+.1*unitrand())*dir(angle)) .. {dir2}pos;
path gs = (0,0){(1,0)} .. cp .. {(0,-1)}cycle;

// draw((0,0)--(3,3));

draw(gs, red);
