#if OPENGL_ON

#pragma once

#include "GL/gli.h"
#include "shape.h"

/* Global variable definitions */

//defined in shape.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

extern bool debug;
extern bool check;

extern const int sample;

int fid = 0;
bool showF=true, showE=true, showV=true, showT=true;
bool all_circumcircle = true;
bool all_vorV = false;

double COM[3], R; //center and radius

#define SQR(x) ((x)*(x))

/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter()   Find the circumcenter of a triangle.                 */
/*                                                                           */
/*  The result is returned both in terms of x-y coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the x-y coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
/*                                                                           */
/*****************************************************************************/

void tricircumcenter(a, b, c, circumcenter, xi, eta)
double a[2];
double b[2];
double c[2];
double circumcenter[2];
double *xi;
double *eta;
{
	double xba, yba, xca, yca;
	double balength, calength;
	double denominator;
	double xcirca, ycirca;

	/* Use coordinates relative to point `a' of the triangle. */
	xba = b[0] - a[0];
	yba = b[1] - a[1];
	xca = c[0] - a[0];
	yca = c[1] - a[1];
	/* Squares of lengths of the edges incident to `a'. */
	balength = xba * xba + yba * yba;
	calength = xca * xca + yca * yca;

	/* Calculate the denominator of the formulae. */
#ifdef EXACT
	/* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
	/*   to ensure a correctly signed (and reasonably accurate) result, */
	/*   avoiding any possibility of division by zero.                  */
	denominator = 0.5 / orient2d(b, c, a);
#else
	/* Take your chances with floating-point roundoff. */
	denominator = 0.5 / (xba * yca - yba * xca);
#endif

	/* Calculate offset (from `a') of circumcenter. */
	xcirca = (yca * balength - yba * calength) * denominator;
	ycirca = (xba * calength - xca * balength) * denominator;
	circumcenter[0] = xcirca;
	circumcenter[1] = ycirca;

	if (xi != (double *)NULL) {
		/* To interpolate a linear function at the circumcenter, define a     */
		/*   coordinate system with a xi-axis directed from `a' to `b' and    */
		/*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
		/*   are computed by Cramer's Rule for solving systems of linear      */
		/*   equations.                                                       */
		*xi = (xcirca * yca - ycirca * xca) * (2.0 * denominator);
		*eta = (ycirca * xba - xcirca * yba) * (2.0 * denominator);
	}
}


void computeCOM_R()
{
    //-------------------------------------------------------------------------
    // compute center of mass and R...
	int total=0;
    int i=0;
    double d=0;
    tVertex  v= vertices;
    COM[0]=COM[1]=COM[2]=0;
    R=0;
    
	do {                                 
	   COM[0]+=v->v[0];
	   COM[1]+=v->v[1];
	   COM[2]+=v->v[2];
	   total++;
	   v = v->next;
	} 
	while ( v != vertices );
   
    for(i=0;i<3;i++) COM[i]/=total;

    v= vertices;
    do {                                 
	   d=SQR(v->v[0]-COM[0])+SQR(v->v[1]-COM[1])+SQR(v->v[2]-COM[2]);
	   if(d>R) R=d;
	   v = v->next;
	} while ( v != vertices );
	
    R=sqrt(R);
}

void DrawTetrahedra()
{
    tTetra   t;
    tVertex  v1, v2, v3, v4;

    if(tetras==NULL) return;

    t = tetras;

    do {

      v1=t->vertex[0];
      v2=t->vertex[1];
      v3=t->vertex[2];
      v4=t->vertex[3];

      //t1
      glVertex3dv(v1->v);
      glVertex3dv(v2->v);
      glVertex3dv(v3->v);
	  
      //t2
      glVertex3dv(v1->v);
      glVertex3dv(v2->v);
      glVertex3dv(v4->v);

      //t3
      glVertex3dv(v2->v);
      glVertex3dv(v3->v);
      glVertex3dv(v4->v);

      //t4
      glVertex3dv(v3->v);
      glVertex3dv(v1->v);
      glVertex3dv(v4->v);

      t = t->next;

    } while ( t != tetras );
}

//draw results using openGL
void DrawTriangles()
{
   
   tFace    f;
   tVertex  v1, v2, v3;
   
   if(faces==NULL) return;

   glEnable( GL_POLYGON_OFFSET_FILL );
   glPolygonOffset( 0.5f, 0.5f );
   glBegin(GL_TRIANGLES);
   f = faces;

   do {           
	
	  v1=f->vertex[0];
	  v2=f->vertex[1];
	  v3=f->vertex[2];

	  glVertex3dv(v1->v);
	  glVertex3dv(v2->v);
      glVertex3dv(v3->v);
      
      f = f->next;
   } while ( f != faces );
   
   glEnd();
   glDisable( GL_POLYGON_OFFSET_FILL );

   glPushAttrib(GL_ENABLE_BIT);
   glColor3f(0,0,0);
   f = faces;

   do {

      v1=f->vertex[0];
      v2=f->vertex[1];
      v3=f->vertex[2];

      glBegin(GL_LINE_LOOP);
      glVertex3dv(v1->v);
      glVertex3dv(v2->v);
      glVertex3dv(v3->v);
      glEnd();

      f = f->next;
   } while ( f != faces );


   glPopAttrib();
}


//draw results using openGL
void DrawTriangles_lifted()
{

	tFace    f;
	tVertex  v1, v2, v3;

	if (faces == NULL) return;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(0.5f, 0.5f);
	glBegin(GL_TRIANGLES);
	f = faces;

	do {

		v1 = f->vertex[0];
		v2 = f->vertex[1];
		v3 = f->vertex[2];

		glVertex3d(v1->v[0], v1->v[1], v1->v[0] * v1->v[0] + v1->v[1] * v1->v[1]);
		glVertex3d(v2->v[0], v2->v[1], v2->v[0] * v2->v[0] + v2->v[1] * v2->v[1]);
		glVertex3d(v3->v[0], v3->v[1], v3->v[0] * v3->v[0] + v3->v[1] * v3->v[1]);

		f = f->next;
	} while (f != faces);

	glEnd();
	glDisable(GL_POLYGON_OFFSET_FILL);

	glPushAttrib(GL_ENABLE_BIT);
	glColor3f(0, 0, 0);
	f = faces;

	do {

		v1 = f->vertex[0];
		v2 = f->vertex[1];
		v3 = f->vertex[2];

		glBegin(GL_LINE_LOOP);
		glVertex3d(v1->v[0], v1->v[1], v1->v[0] * v1->v[0] + v1->v[1] * v1->v[1]);
		glVertex3d(v2->v[0], v2->v[1], v2->v[0] * v2->v[0] + v2->v[1] * v2->v[1]);
		glVertex3d(v3->v[0], v3->v[1], v3->v[0] * v3->v[0] + v3->v[1] * v3->v[1]);
		glEnd();

		f = f->next;
	} while (f != faces);


	glPopAttrib();
}


void DrawEdges()
{   
   tEdge    e;
   tVertex  v1, v2;
   
   if(edges==NULL) return;

   glBegin(GL_LINES);
   /* Edges. */	
   e = edges;
   do {
	  v1=e->endpts[0];
	  v2=e->endpts[1];
	  glVertex3dv(v1->v);
	  glVertex3dv(v2->v);
      e = e->next;
   } while ( e != edges );
   glEnd();
}

void DrawEdges_lifted()
{
	tEdge    e;
	tVertex  v1, v2;

	if (edges == NULL) return;

	glBegin(GL_LINES);
	/* Edges. */
	e = edges;
	do {
		v1 = e->endpts[0];
		v2 = e->endpts[1];
		//glVertex3dv(v1->v);
		//glVertex3dv(v2->v);
		glVertex3d(v1->v[0], v1->v[1], v1->v[0] * v1->v[0] + v1->v[1] * v1->v[1]);
		glVertex3d(v2->v[0], v2->v[1], v2->v[0] * v2->v[0] + v2->v[1] * v2->v[1]);
		e = e->next;
	} while (e != edges);
	glEnd();
}



void DrawVertices(tVertex vlist)
{   
   tVertex  v;
   v = vlist;
   
   if(vlist==NULL) return;

   glBegin(GL_POINTS);
   do {                                 
	  glVertex3dv(v->v);
      v = v->next;
   } while ( v != vlist );
   glEnd();
}

void DrawVertices_lifted(tVertex vlist)
{
	tVertex  v;
	v = vlist;

	if (vlist == NULL) return;

	glBegin(GL_POINTS);
	do {
		glVertex3d(v->v[0], v->v[1], v->v[0] * v->v[0] + v->v[1] * v->v[1]);
		v = v->next;
	} while (v != vlist);
	glEnd();
}

void drawVorVertx()
{
	tFace    f;
	tVertex  v1, v2, v3;

	if (faces == NULL) return;
	f = faces;

	int count = 0;

	glBegin(GL_POINTS);
	do 
	{
		v1 = f->vertex[0];
		v2 = f->vertex[1];
		v3 = f->vertex[2];

		double a[2] = { v1->v[0], v1->v[1] };
		double b[2] = { v2->v[0], v2->v[1] };
		double c[2] = { v3->v[0], v3->v[1] };
		double circumcenter[2];
		tricircumcenter(a, b, c, circumcenter, NULL, NULL);

		//draw circle center
		glVertex3d(circumcenter[0] + a[0], circumcenter[1] + a[1], 0);
		
		f = f->next;
		count++;
	} 
	while (f != faces);
	glEnd();
}

void drawCircumcircle(int fid)
{

	tFace    f;
	tVertex  v1, v2, v3;

	if (faces == NULL) return;

	f = faces;

	int count = 0;

	do {

		if (count == fid || fid<0)
		{

			v1 = f->vertex[0];
			v2 = f->vertex[1];
			v3 = f->vertex[2];

			double a[2] = { v1->v[0], v1->v[1] };
			double b[2] = { v2->v[0], v2->v[1] };
			double c[2] = { v3->v[0], v3->v[1] };
			double circumcenter[2];
			tricircumcenter(a, b, c, circumcenter, NULL, NULL);

			//draw circle center only for random vertex, don't do this when the program is drawing all circles
			if (fid >= 0)
			{
				glBegin(GL_POINTS);
				glVertex3d(circumcenter[0] + a[0], circumcenter[1] + a[1], 0);
				glEnd();
			}

			glBegin(GL_LINE_LOOP);
			double r = sqrt(circumcenter[0] * circumcenter[0] + circumcenter[1] * circumcenter[1]);
			int cir_res = 40;
			const double PI = 3.14159265359;
			double cir_delta = PI * 2 / cir_res;
			for (int i = 0; i < cir_res; i++)
			{
				double x = cos(cir_delta*i)*r;
				double y = sin(cir_delta*i)*r;
				glVertex3d(circumcenter[0] + a[0] + x, circumcenter[1] + a[1] + y, 0);
			}
			glEnd();

			if (fid >= 0) break;
		}

		f = f->next;
		count++;
	} while (f != faces);
	
}

void drawCircumcircle_lifted(int fid)
{

	tFace    f;
	tVertex  v1, v2, v3;

	if (faces == NULL) return;

	f = faces;
	
	int count = 0;

	do {

		if (count == fid || fid<0)
		{

			v1 = f->vertex[0];
			v2 = f->vertex[1];
			v3 = f->vertex[2];

			double a[2] = { v1->v[0], v1->v[1] };
			double b[2] = { v2->v[0], v2->v[1] };
			double c[2] = { v3->v[0], v3->v[1] };
			double circumcenter[2];
			tricircumcenter(a, b, c, circumcenter, NULL, NULL);

			//draw circle center only for random vertex, don't do this when the program is drawing all circles
			if (fid >= 0)
			{
				glBegin(GL_POINTS);
				double cx = circumcenter[0] + a[0];
				double cy = circumcenter[1] + a[1];
				glVertex3d(cx, cy, cy*cy + cx*cx);
				glEnd();
			}

			glBegin(GL_LINE_LOOP);
			double r = sqrt(circumcenter[0] * circumcenter[0] + circumcenter[1] * circumcenter[1]);
			int cir_res = 40;
			const double PI = 3.14159265359;
			double cir_delta = PI * 2 / cir_res;
			for (int i = 0; i < cir_res; i++)
			{
				double x = cos(cir_delta*i)*r + circumcenter[0] + a[0];
				double y = sin(cir_delta*i)*r + circumcenter[1] + a[1];
				glVertex3d(x, y, x*x + y*y);
			}
			glEnd();

			if (fid>=0) break;
		}

		f = f->next;
		count++;

	} 
	while (f != faces);

}

//copied from meshlab
void DisplayBackground(void)
{
	float topcolor[]={1,1,1};
	float bottomcolor[]={0,0,0};
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1,1,-1,1,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLE_STRIP);
			glColor3fv(topcolor);      glVertex2f(-1, 1);
			glColor3fv(bottomcolor);   glVertex2f(-1,-1);
			glColor3fv(topcolor);      glVertex2f( 1, 1);
			glColor3fv(bottomcolor);   glVertex2f( 1,-1);
	glEnd();
	
	glPopAttrib();
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}


void Display(void)
{
    static GLfloat light_position1[] = {  100, 100, 100.0f, 1.0f };
    static GLfloat light_position2[] = { -100, -100, 50.0f, 1.0f };

	//Init Draw
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glMatrixMode(GL_MODELVIEW);
    DisplayBackground();
    
    glPushMatrix();
    glLoadIdentity();

    glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
    glPopMatrix();

	//draw
	glTranslated(-COM[0],-COM[1],-COM[2]);
	
	if (showF)
	{
		glColor3f(1, 1, 1);
		DrawTriangles();
		glColor3f(0.7, 0.7, 1);
		DrawTriangles_lifted();
	}

	if (all_circumcircle)
	{
		glPointSize(5);
		glColor3f(0, 0.5, 0);
		drawCircumcircle(-1);
		glColor3f(0.5, 0.5, 1);
		drawCircumcircle_lifted(-1);
	}
	else{
		glPointSize(5);
		glColor3f(0, 0.5, 0);
		drawCircumcircle(fid);
		glColor3f(0.5, 0.5, 1);
		drawCircumcircle_lifted(fid);
	}
	
	if(showE){
		glColor3f(0,0,0);
		DrawEdges();
		DrawEdges_lifted();
	}
	
	if(showV){
		glPointSize(5);
		glColor3f(1,0,0);
		DrawVertices(vertices);
		glColor3f(0, 0, 1);
		DrawVertices_lifted(vertices);
	}

	if(showT){
	    glEnable( GL_POLYGON_OFFSET_FILL );
	    glPolygonOffset( 0.5f, 0.5f );
	    glBegin(GL_TRIANGLES);
	    glColor3f(1,1,1);
	    DrawTetrahedra();
	    glEnd();
	    glDisable( GL_POLYGON_OFFSET_FILL );

	    glBegin(GL_LINES);
	    glColor3f(0,0,0);
	    DrawTetrahedra();
	    glEnd();
	}

	if (all_vorV)
	{
		glColor3f(0.5, 0.5, 0);
		drawVorVertx();
	}
}

//-----------------------------------------------------------------------------
// other regular openGL callback functions
bool InitGL()
{
    GLfloat Diffuse[] =  { 0.9f, 0.9f, 0.9f, 1.0f };
    GLfloat WhiteLight[] =  { 1.0f, 1.0f, 1.0f, 1.0f };

	// transparent
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    // others
    glEnable( GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

    glClearColor( 1,1,1,0 );

    //Let's have light!

    glMaterialfv(GL_FRONT, GL_DIFFUSE, Diffuse);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);


    glLightfv(GL_LIGHT0,GL_DIFFUSE,WhiteLight);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,WhiteLight);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    return true;
}

void Reshape( int w, int h)
{
    if(w>h)
        glViewport( 0, 0, (GLsizei)w, (GLsizei)w );
    else
        glViewport( 0, 0, (GLsizei)h, (GLsizei)h );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 60, 1, 0.1, 1000 );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}

// show keys for controlling the rendering
void printGUIKeys()
{
	printf("GUI keys:\n");
	printf("r:   select a random circumcircle\n");
	printf("a:   show/hide all circumcircles\n");
	printf("v:   show Voronoi verices\n");

    printf("?:   show this message\n");
    printf("esc: quit\n");
}

//keyboard event function
void Keyboard( unsigned char key, int x, int y )
{	
    // find closest colorPt3D if ctrl is pressed...
    switch( key ){
        case 27: exit(0);
		case 'r': fid = ((rand()*1.0 / RAND_MAX)) * sample; break;
		case 'a': all_circumcircle = !all_circumcircle; break;
		case 'v': all_vorV = !all_vorV; break;

        case '?': printGUIKeys(); break;
    }

    glutPostRedisplay();
}


#endif //OPENGL_ON


void Draw3D(int argc, char ** argv)
{

#if OPENGL_ON

	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH );
	glutInitWindowSize( 800, 800);
	glutInitWindowPosition( 50, 50 );
	glutCreateWindow( "shape" );

	InitGL();
	gliInit();
	gliDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Keyboard);

	/////////////////////////////////////////////////////////////
	computeCOM_R();
	//set camera position
	setCameraPosZ(R*2);
	printGUIKeys();
	/////////////////////////////////////////////////////////////
	gliMainLoop();
	
#else

	fprintf(stderr,"! ERROR: OpenGL is not supported. Please recompile with OPENGL_ON=1 in Makefile\n");
	
#endif //OPENGL_ON

}
