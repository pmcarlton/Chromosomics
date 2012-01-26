//Synaptonemal complex modeling in processing + fvlib
// 020120123_13:47_pmc
// starting new, after canonicalizing scmod7. this time, try nuclear confinement…
//Licensed under CC attribution-sharealike http://creativecommons.org/licenses/by-sa/3.0/

import peasy.org.apache.commons.math.*;
import peasy.*;
import peasy.org.apache.commons.math.geometry.*;
import volatileprototypes.fvlib.*;
// Requires pmcarlton’s branch of fvlib with BehaviorSphereConstraint
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
BehaviorConstantDistance bcd;
BehaviorSpringRelaxation rs;
BehaviorSphereConstraint bsph;
PeasyCam cam;

float Kd = 1.0;
int lstmp=0;
int ww=800;
int hh=600;
float frx=.999;
float ldist=5.0;
float sdist=5.0;
//float stLC=0.5;    //stiffness between Lateral and Central elements
//float stCEL=0.1; //stiffness for consecutive CEs
//float stLEOP=0.5; //stiffness for opposite LEs
float stSKIP=0.241;  //stiffness for every other consecutive LE
float stLAT=0.133;   //stiffness for consecutive LEs
int steps=200;
int len=100;
int it=0;
float msec=0;
PVector origin;
float theRadius=21.0f;
float sphForce=0.00001;
PFont font;
float strokew=3;
boolean sepon=false;
float theta=0.0;
int flip=1;
int pindex=0;
int lindex=0;

void setup() {
  String  _timestamp = "" + System.currentTimeMillis();
  size(ww, hh, OPENGL);
  hint(ENABLE_OPENGL_4X_SMOOTH);  
  hint(ENABLE_NATIVE_FONTS) ;
  //  hint(ENABLE_DEPTH_SORT);
  cam = new PeasyCam(this, 0, 0, 0, 300);       // Camera Setup
  cam.setMinimumDistance(10);
  cam.setMaximumDistance(1500);
  //Make the SC.

  origin=new PVector(0.0f, 0.0f, 0.0f);
  LEs=new Point[2][len];
  allpts=new Point[len*2];
  LinkConsec=new Link[2][len-1];
  LinkEveryOther=new Link[2][len-2];
  //alllinks=new Link[(3*(len))+(3*(len-1))+(2*(len-2))];
  alllinks=new Link[(2*(len-1))+(2*(len-2))];
  for (int i=0;i<len;i++) {
    //     LEs[0][i] = new Point(-(ww/4)+i*ldist,sdist,random(1)/100000,1,0.0,ldist*3);//put the LEs in their own array
    LEs[0][i] = new Point(random(1), random(1), random(1), 1, 0.0, ldist*3);//put the LEs in their own array
    allpts[pindex++]=LEs[0][i];
    //     LEs[1][i] = new Point(-(ww/4)+i*ldist,-sdist,random(1)/100000,1,0.0,ldist*3);
    LEs[1][i] = new Point(random(1), -random(1), random(1), 1, 0.0, ldist*3);
    allpts[pindex++]=LEs[1][i];
  }


  for (int i=0;i<len;i++) {
    if (i>0) {
      LinkConsec[0][i-1]=new Link(LEs[0][i], LEs[0][i-1], stLAT);
      LinkConsec[0][i-1].setC(ldist);
      alllinks[lindex++]=LinkConsec[0][i-1];
      LinkConsec[1][i-1]=new Link(LEs[1][i], LEs[1][i-1], stLAT);
      LinkConsec[1][i-1].setC(ldist);
      alllinks[lindex++]=LinkConsec[1][i-1];
      if (i>1) {
        LinkEveryOther[0][i-2]=new Link(LEs[0][i], LEs[0][i-2], stSKIP);
        LinkEveryOther[0][i-2].setC(ldist*2.0);
        alllinks[lindex++]=LinkEveryOther[0][i-2];
        LinkEveryOther[1][i-2]=new Link(LEs[1][i], LEs[1][i-2], stSKIP);
        LinkEveryOther[1][i-2].setC(ldist*2.0);
        alllinks[lindex++]=LinkEveryOther[1][i-2];
      }
    }
  }


  vi=new IntegratorVerlet(allpts).setF(frx); //initvi
  rs=new BehaviorSpringRelaxation(alllinks);
  bsph=new BehaviorSphereConstraint(allpts, origin, theRadius);
  bsph.setFMult(sphForce);
  font = createFont("Commodore 64", 16);
}

void draw() {

  background(50, 100);
theRadius=30+sin(2*PI*(it%steps)/steps)*10.0;
bsph.setRadius(theRadius);
  for (int j=steps;j>=0;j--) { 
    it++;

    rs.step();
    bsph.step();
    vi.step();
  }

  fill(128, 60);
  stroke(0);
  strokeWeight(0);

  //for (int i=0;i<2;i++) {
  //for (Link l:LinkConsec[i]) {
  //line(l.p1.x,l.p1.y,l.p1.z,l.p2.x,l.p2.y,l.p2.z);
  //}
  //for (Link l:LinkEveryOther[i]){
  //  line(l.p1.x,l.p1.y,l.p1.z,l.p2.x,l.p2.y,l.p2.z);
  //}
  //}
  beginShape(TRIANGLES);  
  for (int i=0;i<2;i++) {
    for (Link l:LinkConsec[i]) {
      vertex(l.p1.x, l.p1.y, l.p1.z);
      vertex(l.p2.x, l.p2.y, l.p2.z);
    }
    //for (Link l:LinkEveryOther[i]){
    //  line(l.p1.x,l.p1.y,l.p1.z,l.p2.x,l.p2.y,l.p2.z);
    //}
  }
  endShape();

  strokeWeight(3);
  beginShape(POINTS);

  for (Point pt : LEs[0]) {
    stroke(255, 0, 120);
    vertex(pt.x, pt.y, pt.z);
  }
  for (Point pt : LEs[1]) {
    stroke(0, 255, 120);
    vertex(pt.x, pt.y, pt.z);
  }
  endShape();


  noFill();
  strokeWeight(0.1);
  stroke(0,200,255,25);
  sphere(theRadius);

  cam.beginHUD();
  textFont(font);
  text("Iteration: " + it, 10, height-10);
  cam.endHUD();
}

void keyPressed() {
  if (keyCode==' ') {
    sepon=!sepon;
  }
}

