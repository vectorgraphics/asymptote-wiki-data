size(0,150);

path ellipse(pair c, real a, real b) {
   return shift(c)*scale(a,b)*unitcircle;
   }

real radius=1.5, eccentricity=2.05;

path E = ellipse((1.2,.2),eccentricity*radius,radius);
path [] P ;
for (int k=0; k<5; ++k) P[k]=rotate(k*72)*E; 

pen [] C = {
   rgb(.45,0,0),rgb(0,.45,0),rgb(0,0,.45),rgb(.30,.30,0),rgb(.30,0,.30),rgb(0,.30,.30)
   };

picture pic;

for (int k=0; k<5; ++k) {
    fill(P[k],C[k]);
    for (int l=0; l<k; ++l) {
       fill(pic, P[k], C[k]+C[l]);
       clip(pic,P[l]);
       add(pic);
       for (int m=0; m<l; ++m) {
          fill(pic, P[k], C[k]+C[l]+C[m]);
          clip(pic,P[l]); clip(pic,P[m]);
          add(pic);
          for (int n=0; n<m; ++n) {
             fill(pic, P[k], C[k]+C[l]+C[m]+C[n]);
             clip(pic,P[l]); clip(pic,P[m]); clip(pic,P[n]);
             add(pic);
          }
       }
    }
}
fill(pic, P[0], C[0]+C[1]+C[2]+C[3]+C[4]);
clip(pic,P[1]); clip(pic,P[2]); clip(pic,P[3]); clip(pic,P[4]);
add(pic);

fill(circle((0,0),.04));

 
