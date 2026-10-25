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

draw((0,0)--(3,3));

draw(randompath(new pair[]{(1,1),(2,1),(2.4,1.8)}, 30));
