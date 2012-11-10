//Synaptonemal complex modeling in processing + fvlib
//020120116_12:26 Peter M. Carlton (pcarlton@icems.kyoto-u.ac.jp)
//Licensed under CC attribution-sharealike http://creativecommons.org/licenses/by-sa/3.0/

/// this one is for getting the coordinates out, for povray/octave image creation
/// used for Biophysical Reviews submission, 2012.11.10


import peasy.org.apache.commons.math.*;
import peasy.*;
import peasy.org.apache.commons.math.geometry.*;
import volatileprototypes.fvlib.*;
import processing.opengl.*;
import processing.video.*;

Point[][] LEs;  //Lateral Elements
Point[] CEs;    //Central Elements
Point[] allpts; //All points, for the Verlet integrator
Link[][] LinkAcross;  //LE-CE and LE-LE links
Link[][] LinkConsec;  //Links between consecutive elements
Link[][] LinkEveryOther;  //Links between every other Lateral (not Central) element
Link[] alllinks;      //All links, for the Spring Relaxation

IntegratorVerlet vi;
BehaviorSpringRelaxation rs;
PeasyCam cam;

float Kd=2100;
float eed=0; //end-to-end distance of whole thing
int lstmp=0;
int ww=800;
int hh=600;
float C=30.0;
float frx=0.99;
float ldist=0.5;
float sdist=0.1;
float stLC=0.5;    //stiffness between Lateral and Central elements
float stCEL=0.05; //stiffness for consecutive CEs
float stLEOP=0.5; //stiffness for opposite LEs
float stSKIP=0.5;  //stiffness for every other consecutive LE
float stLAT=0.5;   //stiffness for consecutive LEs
int steps=200;
int shrinkmod=20000;
int len=51;
int it=0;
float msec=0;
PFont font;
float strokew=3;
boolean sepon=false;
float theta=0.0;
  int flip=1;
PrintWriter output;

void setup() {
int pindex=0;
int lindex=0;
output = createWriter("/tmp/positions.txt"); 
String  _timestamp = "" + System.currentTimeMillis();
  size(ww,hh,OPENGL);
//    size(ww,hh,P3D);
 hint(ENABLE_OPENGL_4X_SMOOTH);  
  cam = new PeasyCam(this, 0,0,0,300);       // Camera Setup
  cam.setMinimumDistance(0);
  cam.setMaximumDistance(1500);
  
  //Make the SC.
LEs=new Point[2][len];
CEs=new Point[len];
allpts=new Point[len*3];
LinkAcross=new Link[3][len];
LinkConsec=new Link[3][len-1];
LinkEveryOther=new Link[2][len-2];
alllinks=new Link[(3*(len))+(3*(len-1))+(2*(len-2))];

   for (int i=0;i<len;i++) {
     LEs[0][i] = new Point(-(2)+i*ldist,sdist,random(1)/100000,1,0.0,ldist*3);//put the LEs in their own array
//     LEs[0][i] = new Point(-(ww/4)+i*ldist,sdist,random(1)/100000,1,0.0,ldist*3);//put the LEs in their own array
     allpts[pindex++]=LEs[0][i];
     LEs[1][i] = new Point(-(2)+i*ldist,-sdist,random(1)/100000,1,0.0,ldist*3);
  //   LEs[1][i] = new Point(-(ww/4)+i*ldist,-sdist,random(1)/100000,1,0.0,ldist*3);
     allpts[pindex++]=LEs[1][i];
     CEs[i]=     new Point(-(2)+i*ldist,0,random(1)/100000,1,0.0,ldist*3); //put the CEs in their own array
// CEs[i]=     new Point(-(ww/4)+i*ldist,0,random(1)/100000,1,0.0,ldist*3); //put the CEs in their own array
    allpts[pindex++]=CEs[i];
   }


  for (int i=0;i<len;i++) {
    LinkAcross[0][i]=new Link(LEs[0][i],CEs[i],stLC);LinkAcross[0][i].setC(sdist);alllinks[lindex++]=LinkAcross[0][i];
    LinkAcross[1][i]=new Link(LEs[0][i],LEs[1][i],stLEOP);LinkAcross[1][i].setC(sdist*2);alllinks[lindex++]=LinkAcross[1][i];
    LinkAcross[2][i]=new Link(LEs[1][i],CEs[i],stLC);LinkAcross[2][i].setC(sdist);alllinks[lindex++]=LinkAcross[2][i];
    if(i>0) {
      LinkConsec[0][i-1]=new Link(LEs[0][i],LEs[0][i-1],stLAT);LinkConsec[0][i-1].setC(ldist);alllinks[lindex++]=LinkConsec[0][i-1];
      LinkConsec[1][i-1]=new Link(CEs[i],CEs[i-1],stCEL);LinkConsec[1][i-1].setC(ldist);alllinks[lindex++]=LinkConsec[1][i-1];
      LinkConsec[2][i-1]=new Link(LEs[1][i],LEs[1][i-1],stLAT);LinkConsec[2][i-1].setC(ldist);alllinks[lindex++]=LinkConsec[2][i-1];
    if(i>1) {
      LinkEveryOther[0][i-2]=new Link(LEs[0][i],LEs[0][i-2],stSKIP);LinkEveryOther[0][i-2].setC(ldist*2.0);alllinks[lindex++]=LinkEveryOther[0][i-2];
      LinkEveryOther[1][i-2]=new Link(LEs[1][i],LEs[1][i-2],stSKIP);LinkEveryOther[1][i-2].setC(ldist*2.0);alllinks[lindex++]=LinkEveryOther[1][i-2];
    }
    }
  }

  vi=new IntegratorVerlet(allpts).setF(frx); //initvi
  rs=new BehaviorSpringRelaxation(alllinks);//.setFast(true);
  font = createFont("Helvetica",16);
}

void draw() {

background(50,100);


  for (int j=steps;j>=0;j--) {
if(0==(it%shrinkmod)) {sepon=(!sepon);}    
    
    it++;
    rs.step();
    vi.step();
    
    if (sepon) {
      eed=sqrt(pow(CEs[0].x-CEs[len-1].x,2) + pow(CEs[0].y-CEs[len-1].y,2) + pow(CEs[0].z-CEs[len-1].z,2));
      for (int i=0;i<(len-1);i++) {
        Link ltmp=LinkConsec[1][i];
if(abs(i-((len-1)/2))==lstmp) 
{

print(it+": "+i+" ");
  ltmp.setC(ltmp.getC()*0.99);
//    ltmp.setC(ltmp.getC()*0.5);
}
    }
    print("\n");
    lstmp++;
    sepon=!sepon;
    }
  }

  fill(255);
  stroke(155);
  strokeWeight(1);
  for (Link l : alllinks) {
    line(C*l.p1.x,C*l.p1.y,C*l.p1.z,C*l.p2.x,C*l.p2.y,C*l.p2.z);
  }

  strokeWeight(7);
  beginShape(POINTS);
  int eee=0;
  for (Point pt : allpts) {
  stroke(127*(eee%3),127*((eee+1)%3),127*((eee+2)%3));
  vertex(C*pt.x,C*pt.y,C*pt.z);
  output.println(C*pt.x+" "+C*pt.y+" "+C*pt.z);
  eee++;
  }
  endShape();
  cam.beginHUD();
  text("Iteration: " + it,10,height-10);
  text("End-to-end distance: " + eed,10,height-30);
//  output.print("eed " + eed + "\n");
  output.flush();
  cam.endHUD();
}

void keyPressed() {
  if (keyCode==' ') {
    sepon=!sepon;
  }
}
