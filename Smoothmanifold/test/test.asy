size(7cm);


srand(0);

path connect (path p, path q)
// Connect `p` and `q` smoothly.
{ return p{dir(p, length(p))}..{dir(q,0)}q; }

pair randomdir (pair dir, real angle)
{ return dir(degrees(dir) + (unitrand()-.5)*angle); }

path randompath (pair[] controlpoints, real angle)
// Like `connect`, but the path has some randomness to it.
{
    if (controlpoints.length < 2) return nullpath;

    pair outdir = randomdir(controlpoints[1]-controlpoints[0], angle);
    guide res = controlpoints[0];

    for (int i = 1; i < controlpoints.length-1; ++i)
    {
        pair indir = randomdir(controlpoints[(i+1)]-controlpoints[i-1], angle);
        res = res{outdir} .. {indir}controlpoints[i];
        outdir = indir;
    }
    
    pair indir = randomdir(controlpoints[(controlpoints.length-1)]-controlpoints[controlpoints.length-2], angle);
    
    return (path)(res{outdir}..{indir}controlpoints[controlpoints.length-1]);
}
// import smoothmanifold;


// element elt = sm1.elements[index1];

pair pos = (1,1);

real angle = 15;
real radius = 1.13;
real range = 30;

// if (dummy(angle))
// { angle = degrees(-elt.labelalign, warn = false); }
// if (dummy(radius))
// { radius = .05*radius(sm1.contour); }

pair dir1 = randomdir(dir(angle-range), range);
pair dir2 = -randomdir(dir(angle+range), range);

path gs = pos{dir1} .. (pos+2*radius*(1+.1*unitrand())*dir(angle)) .. {dir2}pos;

draw((0,0)--(3,3));

draw(gs, red);
// draw(randompath(new pair[]{(1,1),(2,1),(2.4,1.8)}, 30));
