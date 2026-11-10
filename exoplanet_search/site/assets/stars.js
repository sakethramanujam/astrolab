(()=>{const c=document.getElementById('stars');if(!c)return;const x=c.getContext('2d');let s=[];
function r(){c.width=innerWidth*devicePixelRatio;c.height=innerHeight*devicePixelRatio;
s=Array.from({length:Math.round(innerWidth*innerHeight/5000)},()=>({x:Math.random()*c.width,y:Math.random()*c.height,
r:Math.random()*1.2*devicePixelRatio+.2,a:Math.random()*.8+.2,p:Math.random()*6.28}))}
function d(t){x.clearRect(0,0,c.width,c.height);for(const q of s){const b=q.a*(.65+.35*Math.sin(t/1400+q.p));
x.fillStyle=`rgba(220,230,255,${b})`;x.beginPath();x.arc(q.x,q.y,q.r,0,6.283);x.fill()}
if(!matchMedia('(prefers-reduced-motion: reduce)').matches)requestAnimationFrame(d)}
addEventListener('resize',r);r();requestAnimationFrame(d)})();
