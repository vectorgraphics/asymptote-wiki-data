import graph3;

settings.outformat="html";
settings.render=4;

size(10cm);
defaultpen(fontsize(5));

currentprojection=perspective(
    camera=(21.6,-21.2,11.8),
    up=(-0.006,0.008,0.038),
    target=(3.5,3.4,2.9)
);

interaction constantsize = settings.autobillboard
    ? interaction(1,true)
    : interaction(0,true);


// ============================================================
// 1. Initial rigid-body position
// ============================================================

triple[] sp = {
    (2,3,3.1),
    (6,3,3.1),
    (6,7,3.1),
    (2,7,3.1)
};

triple T1c = (sp[0]+sp[1]+sp[2]+sp[3])/4;


// ============================================================
// 2. Static background: ground plane
// ============================================================

real groundMinX = -2;
real groundMaxX = 10;
real groundMinY = -1;
real groundMaxY = 11;
real groundZ = 0;

triple G0 = (groundMinX,groundMinY,groundZ);
triple G1 = (groundMaxX,groundMinY,groundZ);
triple G2 = (groundMaxX,groundMaxY,groundZ);
triple G3 = (groundMinX,groundMaxY,groundZ);

surface ground = surface(G0--G1--G2--G3--cycle);
draw(ground,lightgray+opacity(0.25));

real gridStep = 1.0;

// X-direction grid lines
for(real x=groundMinX; x<=groundMaxX+1e-6; x+=gridStep) {
    draw((x,groundMinY,groundZ+0.01)--(x,groundMaxY,groundZ+0.01),
         gray(0.75)+0.35);
}

// Y-direction grid lines
for(real y=groundMinY; y<=groundMaxY+1e-6; y+=gridStep) {
    draw((groundMinX,y,groundZ+0.01)--(groundMaxX,y,groundZ+0.01),
         gray(0.75)+0.35);
}


// ============================================================
// 3. Front-end algorithm and UI control logic
// ============================================================

javascript("
window.rigidParams={
    tau:0.0,
    tx:2.0,ty:0.5,tz:1.0,
    rx:5.0,ry:10.0,rz:0.0
};

window.T1_CENTER=[4.0,5.0,3.1];


// ------------------------------------------------------------
// 1. Standard RPY composite rotation matrix
// ------------------------------------------------------------

function getTargetRotationMatrix() {
    let rx=window.rigidParams.rx*Math.PI/180;
    let ry=window.rigidParams.ry*Math.PI/180;
    let rz=window.rigidParams.rz*Math.PI/180;

    let cphi=Math.cos(rx);
    let sphi=Math.sin(rx);
    let cthe=Math.cos(ry);
    let sthe=Math.sin(ry);
    let cpsi=Math.cos(rz);
    let spsi=Math.sin(rz);

    return [
        cpsi*cthe,
        spsi*cthe,
        -sthe,

        cpsi*sthe*sphi-spsi*cphi,
        spsi*sthe*sphi+cpsi*cphi,
        cthe*sphi,

        cpsi*sthe*cphi+spsi*sphi,
        spsi*sthe*cphi-cpsi*sphi,
        cthe*cphi
    ];
}


// ------------------------------------------------------------
// 2. Generate a 3x3 rotation matrix using Rodrigues' formula
// ------------------------------------------------------------

function getAxisAngleMatrix(k,Phi) {
    let c=Math.cos(Phi);
    let s=Math.sin(Phi);
    let v=1-c;

    let kx=k[0];
    let ky=k[1];
    let kz=k[2];

    return [
        kx*kx*v+c,
        kx*ky*v+kz*s,
        kx*kz*v-ky*s,

        kx*ky*v-kz*s,
        ky*ky*v+c,
        ky*kz*v+kx*s,

        kx*kz*v+ky*s,
        ky*kz*v-kx*s,
        kz*kz*v+c
    ];
}


// ------------------------------------------------------------
// 3. Core rigid-body transformation
// ------------------------------------------------------------

window.J_Rigid=function(p,t) {

    if(!p || p.length<3)
        return [0,0,0];

    let R_end=getTargetRotationMatrix();

    let trace=R_end[0]+R_end[4]+R_end[8];

    let Phi_total=Math.acos(
        Math.min(Math.max((trace-1)/2,-1),1)
    );

    let k_axis=[0,0,1];

    if(Math.abs(Phi_total)>1e-6) {
        let sinPhi2=2*Math.sin(Phi_total);

        k_axis=[
            (R_end[5]-R_end[7])/sinPhi2,
            (R_end[6]-R_end[2])/sinPhi2,
            (R_end[1]-R_end[3])/sinPhi2
        ];
    }

    let tau=window.rigidParams.tau;
    let cur_Phi=tau*Phi_total;
    let R_mat=getAxisAngleMatrix(k_axis,cur_Phi);

    // Rotate around the rigid-body center
    let px=p[0]-window.T1_CENTER[0];
    let py=p[1]-window.T1_CENTER[1];
    let pz=p[2]-window.T1_CENTER[2];

    let rx=R_mat[0]*px+R_mat[3]*py+R_mat[6]*pz;
    let ry=R_mat[1]*px+R_mat[4]*py+R_mat[7]*pz;
    let rz=R_mat[2]*px+R_mat[5]*py+R_mat[8]*pz;

    return [
        rx+tau*window.rigidParams.tx+window.T1_CENTER[0],
        ry+tau*window.rigidParams.ty+window.T1_CENTER[1],
        rz+tau*window.rigidParams.tz+window.T1_CENTER[2]
    ];
};


// ============================================================
// 4. User interface
// ============================================================

window.addEventListener('load',function() {

    let hide=document.createElement('style');

    hide.textContent=
        '.asy-player-controls,#asy-slider {'+
        'display:none !important;}';

    document.head.appendChild(hide);

    let panel=document.createElement('div');

    panel.style.cssText=
        'position:fixed;'+
        'top:15px;'+
        'left:15px;'+
        'z-index:9999;'+
        'background:rgba(255,255,255,0.92);'+
        'padding:12px;'+
        'border-radius:8px;'+
        'font-family:sans-serif;'+
        'box-shadow:0 4px 12px rgba(0,0,0,0.25);'+
        'width:280px;'+
        'font-size:12px;';

    panel.innerHTML=
        '<b style=\"display:block;'+
        'margin-bottom:8px;'+
        'border-bottom:1px solid #ccc;'+
        'padding-bottom:4px;'+
        'font-size:13px;\">'+
        '6-DoF Rigid-Body Transformation</b>'+

        '<div style=\"display:flex;'+
        'justify-content:space-between;\">'+
        '<label>Progress (tau):</label>'+
        '<span id=\"val-tau\">0.00</span>'+
        '</div>'+

        '<input type=\"range\" id=\"tauSlider\" '+
        'style=\"width:100%;margin:2px 0 8px 0;\" '+
        'min=\"0\" max=\"1\" step=\"0.01\" value=\"0\">'+

        '<div style=\"font-weight:bold;'+
        'margin-top:4px;color:#333;\">'+
        'Translation:</div>'+

        '<div style=\"display:grid;'+
        'grid-template-columns:1fr 1fr 1fr;'+
        'gap:4px;margin-top:2px;\">'+

        '<div>X: <input type=\"number\" '+
        'id=\"inp-tx\" value=\"2.0\" step=\"0.1\" '+
        'style=\"width:45px;\"></div>'+

        '<div>Y: <input type=\"number\" '+
        'id=\"inp-ty\" value=\"0.5\" step=\"0.1\" '+
        'style=\"width:45px;\"></div>'+

        '<div>Z: <input type=\"number\" '+
        'id=\"inp-tz\" value=\"1.0\" step=\"0.1\" '+
        'style=\"width:45px;\"></div>'+

        '</div>'+

        '<div style=\"font-weight:bold;'+
        'margin-top:8px;color:#333;\">'+
        'RPY Orientation Angles (degrees):</div>'+

        '<div style=\"display:grid;'+
        'grid-template-columns:1fr 1fr 1fr;'+
        'gap:4px;margin-top:2px;\">'+

        '<div>Rx: <input type=\"number\" '+
        'id=\"inp-rx\" value=\"5.0\" step=\"1\" '+
        'style=\"width:45px;\"></div>'+

        '<div>Ry: <input type=\"number\" '+
        'id=\"inp-ry\" value=\"10.0\" step=\"1\" '+
        'style=\"width:45px;\"></div>'+

        '<div>Rz: <input type=\"number\" '+
        'id=\"inp-rz\" value=\"0.0\" step=\"1\" '+
        'style=\"width:45px;\"></div>'+

        '</div>'+

        '<button id=\"playBtn\" '+
        'style=\"width:100%;margin-top:10px;'+
        'padding:5px;cursor:pointer;'+
        'font-weight:bold;\">'+
        'Play / Pause</button>';

    document.body.appendChild(panel);

    let slider=document.getElementById('tauSlider');
    let valSpan=document.getElementById('val-tau');
    let playBtn=document.getElementById('playBtn');
    let isPlaying=false;

    let refreshScene=function() {
        if(window.updateScene)
            window.updateScene();
        else if(window.asy && window.asy.render)
            window.asy.render();
    };

    let bindInput=function(id,key) {
        let el=document.getElementById(id);

        el.oninput=el.onchange=function() {
            window.rigidParams[key]=parseFloat(this.value)||0;
            refreshScene();
        };
    };

    bindInput('inp-tx','tx');
    bindInput('inp-ty','ty');
    bindInput('inp-tz','tz');
    bindInput('inp-rx','rx');
    bindInput('inp-ry','ry');
    bindInput('inp-rz','rz');

    slider.oninput=function() {
        window.rigidParams.tau=parseFloat(this.value);
        valSpan.innerText=parseFloat(this.value).toFixed(2);
        refreshScene();
    };

    playBtn.onclick=function() {
        isPlaying=!isPlaying;
        if(isPlaying)
            animate();
    };

    function animate() {
        if(!isPlaying)
            return;

        let val=parseFloat(slider.value)+0.01;

        if(val>1.0)
            val=0.0;

        slider.value=val;
        window.rigidParams.tau=val;
        valSpan.innerText=val.toFixed(2);

        refreshScene();
        requestAnimationFrame(animate);
    }
});
");


// ============================================================
// 5. Define the moving rigid body: rectangular cuboid
// ============================================================

real boxX=4.0;
real boxY=2.5;
real boxZ=1.2;

triple C=T1c;

triple[] box={
    // Bottom vertices
    C+(-boxX/2,-boxY/2,-boxZ/2),
    C+( boxX/2,-boxY/2,-boxZ/2),
    C+( boxX/2, boxY/2,-boxZ/2),
    C+(-boxX/2, boxY/2,-boxZ/2),

    // Top vertices
    C+(-boxX/2,-boxY/2, boxZ/2),
    C+( boxX/2,-boxY/2, boxZ/2),
    C+( boxX/2, boxY/2, boxZ/2),
    C+(-boxX/2, boxY/2, boxZ/2)
};


// ============================================================
// 6. Dynamic rigid body
// ============================================================

beginTransform("function(x,t){ return J_Rigid(x,t); }",10);

    // Bottom edges
    draw(box[0]--box[1],magenta+1.5);
    draw(box[1]--box[2],magenta+1.5);
    draw(box[2]--box[3],magenta+1.5);
    draw(box[3]--box[0],magenta+1.5);

    // Top edges
    draw(box[4]--box[5],magenta+1.5);
    draw(box[5]--box[6],magenta+1.5);
    draw(box[6]--box[7],magenta+1.5);
    draw(box[7]--box[4],magenta+1.5);

    // Vertical edges
    draw(box[0]--box[4],magenta+1.5);
    draw(box[1]--box[5],magenta+1.5);
    draw(box[2]--box[6],magenta+1.5);
    draw(box[3]--box[7],magenta+1.5);

    real axLen=1.5;

    triple O=T1c;
    triple X=T1c+(axLen,0,0);
    triple Y=T1c+(0,axLen,0);
    triple Z=T1c+(0,0,axLen);

    // Coordinate axes
    draw(O--X,red+1.2);
    draw(O--Y,green+1.2);
    draw(O--Z,blue+1.2);

    // Coordinate labels
    label("$O$",O,SW);
    label("$X$",X,E);
    label("$Y$",Y,N);
    label("$Z$",Z,N);

endTransform();


// ============================================================
// 7. Target pose
// ============================================================

beginTransform(
    "function(x,t){ "+
    "let sav=window.rigidParams.tau;"+
    "window.rigidParams.tau=1.0;"+
    "let res=J_Rigid(x,t);"+
    "window.rigidParams.tau=sav;"+
    "return res; }",
    10
);

    // Bottom edges
    draw(box[0]--box[1],dashed+orange+1.2);
    draw(box[1]--box[2],dashed+orange+1.2);
    draw(box[2]--box[3],dashed+orange+1.2);
    draw(box[3]--box[0],dashed+orange+1.2);

    // Top edges
    draw(box[4]--box[5],dashed+orange+1.2);
    draw(box[5]--box[6],dashed+orange+1.2);
    draw(box[6]--box[7],dashed+orange+1.2);
    draw(box[7]--box[4],dashed+orange+1.2);

    // Vertical edges
    draw(box[0]--box[4],dashed+orange+1.2);
    draw(box[1]--box[5],dashed+orange+1.2);
    draw(box[2]--box[6],dashed+orange+1.2);
    draw(box[3]--box[7],dashed+orange+1.2);

endTransform();