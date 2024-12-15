import React, { useState, useRef, useEffect } from 'react';
import { createRoot } from 'react-dom/client';

function solveHomography(p0, p1, p2, p3) {
  const U = [0,1,1,0];
  const V = [0,0,1,1];
  const X = [p0.x, p1.x, p2.x, p3.x];
  const Y = [p0.y, p1.y, p2.y, p3.y];

  const A = [];
  const B = [];

  for (let i=0; i<4; i++){
    const u = U[i], v=V[i], x=X[i], y=Y[i];
    A.push([u, v, 1, 0, 0, 0, -u*x, -v*x]);
    B.push(x);
    A.push([0, 0, 0, u, v, 1, -u*y, -v*y]);
    B.push(y);
  }

  function gaussianElimination(m, b) {
    const n = m.length;
    for (let i=0; i<n; i++){
      let pivot = i;
      for (let r=i+1; r<n; r++){
        if (Math.abs(m[r][i]) > Math.abs(m[pivot][i])) pivot=r;
      }
      if (pivot!==i) {
        [m[i], m[pivot]]=[m[pivot], m[i]];
        [b[i], b[pivot]]=[b[pivot], b[i]];
      }
      const pivotVal=m[i][i];
      if (Math.abs(pivotVal)<1e-14) continue;
      for (let r=i+1; r<n; r++){
        const f=m[r][i]/pivotVal;
        b[r]-=f*b[i];
        for (let c=i;c<n;c++){
          m[r][c]-=f*m[i][c];
        }
      }
    }
    const x=new Array(n).fill(0);
    for (let i=n-1;i>=0;i--){
      let sum=b[i];
      for (let c=i+1;c<n;c++){
        sum-=m[i][c]*x[c];
      }
      x[i]=m[i][i]===0?0:sum/m[i][i];
    }
    return x;
  }

  const h = gaussianElimination(A.map(r=>r.slice()), B.slice());

  return {
    h11: h[0], h12: h[1], h13: h[2],
    h21: h[3], h22: h[4], h23: h[5],
    h31: h[6], h32: h[7], h33: 1
  };
}

function invertHomography(H) {
  const {h11,h12,h13,h21,h22,h23,h31,h32,h33}=H;
  const m = [
    [h11,h12,h13],
    [h21,h22,h23],
    [h31,h32,h33]
  ];
  const det = 
    m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) -
    m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +
    m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);

  if (Math.abs(det)<1e-14) throw new Error("Non-invertible homography");

  const inv = [
    [
      (m[1][1]*m[2][2]-m[1][2]*m[2][1])/det,
      (m[0][2]*m[2][1]-m[0][1]*m[2][2])/det,
      (m[0][1]*m[1][2]-m[0][2]*m[1][1])/det
    ],
    [
      (m[1][2]*m[2][0]-m[1][0]*m[2][2])/det,
      (m[0][0]*m[2][2]-m[0][2]*m[2][0])/det,
      (m[0][2]*m[1][0]-m[0][0]*m[1][2])/det
    ],
    [
      (m[1][0]*m[2][1]-m[1][1]*m[2][0])/det,
      (m[0][1]*m[2][0]-m[0][0]*m[2][1])/det,
      (m[0][0]*m[1][1]-m[0][1]*m[1][0])/det
    ]
  ];

  return {
    h11: inv[0][0], h12: inv[0][1], h13: inv[0][2],
    h21: inv[1][0], h22: inv[1][1], h23: inv[1][2],
    h31: inv[2][0], h32: inv[2][1], h33: inv[2][2]
  };
}

function matmul3(A,B) {
  const R = [];
  for (let i=0;i<3;i++){
    R[i]=[];
    for (let j=0;j<3;j++){
      R[i][j]=A[i][0]*B[0][j]+A[i][1]*B[1][j]+A[i][2]*B[2][j];
    }
  }
  return R;
}

// Circle conic in u,v
function circleConic() {
  return [
    [1,    0,   -0.5],
    [0,    1,   -0.5],
    [-0.5,-0.5, 0.25]
  ];
}

function transformConic(C_u, H) {
  const Hinv = invertHomography(H);
  const HinvMat = [
    [Hinv.h11,Hinv.h12,Hinv.h13],
    [Hinv.h21,Hinv.h22,Hinv.h23],
    [Hinv.h31,Hinv.h32,Hinv.h33]
  ];
  const HinvT = [
    [Hinv.h11,Hinv.h21,Hinv.h31],
    [Hinv.h12,Hinv.h22,Hinv.h32],
    [Hinv.h13,Hinv.h23,Hinv.h33]
  ];

  const temp = matmul3(HinvT, C_u);
  const C_x = matmul3(temp, HinvMat);
  return C_x;
}

function ellipseParams(C) {
  const a = C[0][0];
  const b = C[1][1];
  const d = C[0][1]; 
  const g = C[0][2];
  const f = C[1][2];
  const c = C[2][2];

  // Solve for center
  const det = a*b - d*d;
  let x0=0,y0=0;
  if (Math.abs(det)>1e-14) {
    x0 = (b*(-g)-d*(-f))/det;
    y0 = (-d*(-g)+a*(-f))/det;
  }

  const cPrime = a*x0*x0+2*d*x0*y0+b*y0*y0+2*g*x0+2*f*y0+c;
  let theta=0;
  if (Math.abs(a-b)>1e-14) {
    theta = 0.5*Math.atan2(2*d,(a - b));
  } else {
    if (Math.abs(d)>1e-14) theta=Math.PI/4; 
    else theta=0;
  }

  const cosT = Math.cos(theta);
  const sinT = Math.sin(theta);
  const A11p = a*cosT*cosT + 2*d*cosT*sinT + b*sinT*sinT;
  const A22p = a*sinT*sinT - 2*d*sinT*cosT + b*cosT*cosT;

  const rx = Math.sqrt(Math.abs(-cPrime/A11p));
  const ry = Math.sqrt(Math.abs(-cPrime/A22p));
  const angleDeg = theta*180/Math.PI;

  return {
    cx: x0,
    cy: y0,
    rx,
    ry,
    angle: angleDeg
  };
}

// Convexity check
function isConvexPolygon(points) {
  const n=points.length;
  if (n<4) return true;
  let sign=0;
  for (let i=0;i<n;i++){
    const p0=points[i];
    const p1=points[(i+1)%n];
    const p2=points[(i+2)%n];
    const dx1=p1.x-p0.x;
    const dy1=p1.y-p0.y;
    const dx2=p2.x-p1.x;
    const dy2=p2.y-p1.y;
    const cross=dx1*dy2 - dy1*dx2;
    // Debug log cross products
    // console.log(`Edge ${i}: cross=${cross}`);
    if (cross!==0){
      const s=cross>0?1:-1;
      if (sign===0) sign=s;
      else if(sign!==s) {
        return false;
      }
    }
  }
  return true;
}

function QuadrilateralEditor() {
  const initialPoints = [
    { x: 100, y: 100 },
    { x: 300, y: 100 },
    { x: 300, y: 300 },
    { x: 100, y: 300 }
  ];

  const [points, setPoints] = useState(initialPoints);
  const [dragIndex, setDragIndex] = useState(null);
  const [offset, setOffset] = useState({x:0,y:0});
  const [convex, setConvex] = useState(true);
  const svgRef = useRef(null);

  const handleMouseDown = (e, i) => {
    e.stopPropagation();
    const rect = svgRef.current.getBoundingClientRect();
    const px = e.clientX - rect.left;
    const py = e.clientY - rect.top;
    setDragIndex(i);
    setOffset({x: px - points[i].x, y: py - points[i].y});
  };

  const handleMouseMove = (e) => {
    if (dragIndex!==null) {
      const rect = svgRef.current.getBoundingClientRect();
      const px = e.clientX - rect.left;
      const py = e.clientY - rect.top;
      const candidatePoints = points.map((p,idx)=>
        idx===dragIndex?{x:px - offset.x, y:py - offset.y}:p
      );
      const c = isConvexPolygon(candidatePoints);
      console.log("Checking convexity for candidate:", candidatePoints, "Result:", c);
      setConvex(c);
      if (c) {
        setPoints(candidatePoints);
      }
      // If not convex, we don't update points, effectively blocking movement.
    }
  };

  const handleMouseUp = () => {
    if (dragIndex!==null) {
      setDragIndex(null);
    }
  };

  useEffect(()=>{
    const handleUp=()=>handleMouseUp();
    const handleMove=(e)=>handleMouseMove(e);
    window.addEventListener('mouseup',handleUp);
    window.addEventListener('mousemove',handleMove);
    return ()=>{
      window.removeEventListener('mouseup',handleUp);
      window.removeEventListener('mousemove',handleMove);
    }
  });

  const [p0,p1,p2,p3]=points;
  const d1=[p0,p2];
  const d2=[p1,p3];

  const m01={x:(p0.x+p1.x)/2,y:(p0.y+p1.y)/2};
  const m12={x:(p1.x+p2.x)/2,y:(p1.y+p2.y)/2};
  const m23={x:(p2.x+p3.x)/2,y:(p2.y+p3.y)/2};
  const m30={x:(p3.x+p0.x)/2,y:(p3.y+p0.y)/2};

  let h, C_u, C_x, ellipseP;
  try {
    h = solveHomography(p0,p1,p2,p3);
    C_u = circleConic();
    C_x = transformConic(C_u, h);
    ellipseP = ellipseParams(C_x);
  } catch(err) {
    console.error("Error computing ellipse:", err);
  }

  return (
    <svg 
      ref={svgRef} 
      width={500} 
      height={500} 
      style={{border:'1px solid black', userSelect:'none'}}
      onMouseDown={(e)=>e.stopPropagation()}>

      {/* Debug text */}
      <text x="10" y="20" fill="black" fontSize="16">
        Convex: {convex ? "YES" : "NO"}
      </text>

      {/* Quadrilateral */}
      <polygon 
        points={points.map(p=>`${p.x},${p.y}`).join(' ')}
        fill="none" 
        stroke="black" 
        strokeWidth="2"/>

      {/* Diagonals */}
      <line x1={d1[0].x} y1={d1[0].y} x2={d1[1].x} y2={d1[1].y} stroke="red" strokeDasharray="4,2"/>
      <line x1={d2[0].x} y1={d2[0].y} x2={d2[1].x} y2={d2[1].y} stroke="red" strokeDasharray="4,2"/>

      {/* Mid-segment lines */}
      <line x1={m01.x} y1={m01.y} x2={m23.x} y2={m23.y} stroke="blue" strokeDasharray="4,4"/>
      <line x1={m12.x} y1={m12.y} x2={m30.x} y2={m30.y} stroke="blue" strokeDasharray="4,4"/>

      {/* Ellipse (if computed) */}
      {ellipseP && !isNaN(ellipseP.rx) && !isNaN(ellipseP.ry) && ellipseP.rx>0 && ellipseP.ry>0 && (
        <ellipse 
          cx={ellipseP.cx} 
          cy={ellipseP.cy} 
          rx={ellipseP.rx} 
          ry={ellipseP.ry} 
          fill="rgba(0,255,0,0.2)" 
          stroke="green"
          transform={`rotate(${ellipseP.angle},${ellipseP.cx},${ellipseP.cy})`}
        />
      )}

      {/* Control points */}
      {points.map((p,i)=>(
        <circle 
          key={i}
          cx={p.x} 
          cy={p.y} 
          r={8}
          fill="orange"
          stroke="black"
          style={{cursor:'grab'}}
          onMouseDown={(e)=>handleMouseDown(e,i)}
        />
      ))}
    </svg>
  );
}

const root = createRoot(document.getElementById('root'));
root.render(<QuadrilateralEditor />);
