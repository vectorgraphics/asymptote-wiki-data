size(7cm);

pair cp = (2.182992367,0.584931041); // path self-intersects
// pair cp = (2.18299236,0.58493104); // path doesn't self-intersect

path p = (0,0){(1,0)} .. cp .. {(0,-1)}cycle;

draw(p);
