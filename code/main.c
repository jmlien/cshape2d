
#include <limits.h>
#include <assert.h>
#include "draw.h"
#include "shape.h"

/*-------------------------------------------------------------------*/
// global data
const int width  = 1;
const int height = 1;
const int sample = 1000;

/*-------------------------------------------------------------------*/
//defined in shape.h/.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

/*-------------------------------------------------------------------*/
#include "qhull.h"
#include "poly.h"
#include "qset.h"

tEdge QHullCreateEdge(tVertex v1[3], tVertex v2[3])
{
    int i, j, in, jn;
    bool found = false;
	tEdge e=NULL;

    for(i=0;i<3;i++)
	{
        in=i+1;
        if(in==3) in=0;
        for(j=0;j<3;j++){
            jn=j+1;
            if(jn==3) jn=0;
            if(v1[i]==v2[j] && v1[in]==v2[jn]){ found=true; break; }
            if(v1[in]==v2[j] && v1[i]==v2[jn]){ found=true; break; }
        }

        if(found) break;
    }

	assert(found);

    e=MakeNullEdge();
    e->endpts[0]=v1[i];
    e->endpts[1]=v1[in];

    return e;
}

void Delaunay()
{
	tVertex  ptr_v;
	tVertex * all_v=NULL;
	int vsize=0;
	int id=0;
	
	//global varibles for qhull
	static char * options=(char *)"delaunay QJ Pp";
    int curlong, totlong;
	coordT * pt=NULL;
    facetT *facet=NULL;
    vertexT *vertex=NULL;
    vertexT **vertexp=NULL;
    facetT *neighbor, **neighborp;
    int vid=0;
	tFace  face;

    //count number of points
    ptr_v=vertices;
	do {                                 
	   vsize++;
	   ptr_v = ptr_v->next;
	} while ( ptr_v != vertices );
    
    //allocate memory
    pt=(coordT*)calloc(vsize*3,sizeof(coordT)); //each point will have three coord
    all_v=(tVertex*)calloc(vsize,sizeof(tVertex));
    assert(pt && all_v);

    //copy points
    ptr_v=vertices;
	do {                                 
	   pt[id++]=ptr_v->v[0];
	   pt[id++]=ptr_v->v[1];
	   pt[id++] = ptr_v->v[0] * ptr_v->v[0] + ptr_v->v[1] * ptr_v->v[1];
	   all_v[ptr_v->vnum]=ptr_v;
	   ptr_v = ptr_v->next;
	} 
	while ( ptr_v != vertices );

    //using qhull

    qh_init_A(stdin, stdout, stderr, 0, NULL);

    qh DELAUNAY= True;     /* 'd'   */
    //qh SCALElast= True;    /* 'Qbb' */
    //qh KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */

    qh_initflags(options);
    qh_init_B (pt, vsize, 3, false);
    qh_qhull();
    qh_check_output();

    //loop through all faces
    FORALLfacets 
	{    
		if (facet->upperdelaunay) continue;

        face = MakeNullFace(); //make a face

        //get vertices of facet
        //loop through each vertex
        vid=0;
        FOREACHvertex_(facet->vertices)
        {
            //get the id of the vertex
            face->vertex[vid++]=all_v[qh_pointid(vertex->point)];
        }

        FOREACHneighbor_(facet)
        {
            if(facet<neighbor){

                tVertex vertices[3];
                vid=0;
                FOREACHvertex_(neighbor->vertices)
                {
                    //get vertex
                    vertices[vid++]=all_v[qh_pointid(vertex->point)];
                }

                QHullCreateEdge(face->vertex, vertices);
            }
        }//FOREACHneighbor_

    }
    
    //not used
    free(pt); 
    free(all_v);
    pt=NULL;
    all_v=NULL;
    
    //free mem
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
}

void CreateVertices()
{
	tVertex  v;
	int	    vnum = 0;

	srand(time(NULL));

	for(int i=0;i<sample;i++)  
	{
		v = MakeNullVertex();
		v->v[X] = ((rand()*1.0 / RAND_MAX) - 0.5) * width;
		v->v[Y] = ((rand()*1.0 / RAND_MAX) - 0.5) * height;
		v->v[Z] = 0;
		v->vnum = vnum++;
		if ((fabs(v->v[X]) > SAFE) || (fabs(v->v[Y]) > SAFE) || (fabs(v->v[Z]) > SAFE)) {
			printf("Coordinate of vertex below might be too large: run with -d flag\n");
			PrintPoint(v);
		}
	}/*end while*/
}


/*-------------------------------------------------------------------*/
int main( int argc, char *argv[] )
{
   
   CreateVertices();

   Delaunay(); //this simply shows an example of how qhull can be used to build CH
      
   if(vertices==NULL && faces==NULL && edges==NULL && tetras==NULL){
       printf("! Error: empty tetras, vertices, faces and edges\n"); 
       return 1;
   }
   
   Draw3D(argc,argv);
   
   return 0;
}

