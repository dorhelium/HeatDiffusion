//260761484
//Doreen He
//Doreen He

package comp557.a3;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

import java.lang.Math;

/**
 * Half edge data structure.
 * Maintains a list of faces (i.e., one half edge of each) to allow
 * for easy display of geometry.
 */
public class HEDS {
	
    /** List of faces */
    ArrayList<Face> faces = new ArrayList<Face>();
    
    /** List of vertices */
    ArrayList<Vertex> vertices;

    /** Convenience member for keeping track of half edges you make or need */
    Map<String,HalfEdge> halfEdges = new TreeMap<String,HalfEdge>();
    
    /**
     * Builds a half edge data structure from the polygon soup   
     * @param soup
     */
    public HEDS( PolygonSoup soup ){
        vertices = soup.vertexList;
        
       
      
        
        
        for ( int[] f : soup.faceList ){        
        	// TODO: 2 Build the half edge data structure from the polygon soup, triangulating non-triangular faces
        	/*
        	 * Create the half edge and face objects to make up your half edge data structure. 
        	 * See the createHalfEdge method in HEDS which will help you create twin pointers 
        	 * in your half edge. Triangulate any faces in the soup which do not have 3 sides. 
        	 * Take care to ensure that you set all the necessary pointers, 
        	 * and test your half edge data structure by walking around the mesh. 
        	 * Set "draw test half edge" in the controls to true, 
        	 * and then use space, n, and w keys on the keyboard. 
        	 * You may need to rotate your mesh or view in wireframe to see the half edge! 
        	 */
        	
        	// triangulate faces
        	ArrayList<int[]> subfaces = new ArrayList<int[]> ();
        	
        	if (f.length>3) {
        		for (int i=1; i<=f.length-2; i++) {
        			int[] subF = {f[0], f[i], f[i+1]};
        			subfaces.add(subF);
        		}
        	}else {
        		subfaces.add(f);
        	}
        	
        	
        	for (int[] face: subfaces) {

	        	try{
	        		//half edge data structure
		        	HalfEdge start = createHalfEdge (soup, face[face.length-1], face[0]);
		        	HalfEdge prev = start;
		        	
		        	HalfEdge currentEdge = null;
		        	for (int i=0; i<=face.length-2; i++) {
			        	currentEdge = createHalfEdge(soup, face[i], face[i+1]);
			        	//vertices.add(currentEdge.head);
			        	prev.next = currentEdge;
			        	prev = currentEdge;  	
		        	}
		        	
		        	prev.next = start;
		        	Face newFace = new Face(start);
		        	faces.add(newFace);
	        	
	        	}catch (Exception e) {
	        		e.printStackTrace();
	        	}        	
        	
        	}
	        
        }
        
        
        
        
        
        
        
        // TODO: 3 Compute vertex normals
        for (Vertex v: vertices) {

        	ArrayList<Vector3d> adjNorms = new ArrayList<Vector3d>();
        	
        	HalfEdge start = v.he;
        	HalfEdge loop = start;
        	do {
        		adjNorms.add(loop.leftFace.n);
        		loop = loop.twin.next.next;

        	}while(!(loop.equals(start)));
  
        	
        	Vector3d aveNorm = new Vector3d();
        	
        	for (Vector3d norm: adjNorms) {
        		aveNorm.x = aveNorm.x + norm.x;
        		aveNorm.y = aveNorm.y + norm.y;
        		aveNorm.z = aveNorm.z + norm.z;
        		
        	}
        	
        	v.n = aveNorm;
  	
        	
        }
        
        
        
        
    }
    
    
    
    
    
    
    
    /**
     * Helper function for creating a half edge, and pairing it up with its twin if the
     * twin half edge has already been created.
     * @param soup 
     * @param i tail vertex index
     * @param j head vertex index
     * @return the half edge, paired with its twin if the twin was created.
     * @throws Exception
     */
    private HalfEdge createHalfEdge( PolygonSoup soup, int i, int j ) throws Exception {
        String p = i+","+j;
        if ( halfEdges.containsKey( p ) ){
            throw new Exception("non orientable manifold");
        }
        String twin = j+","+i;
        HalfEdge he = new HalfEdge();
        he.head = soup.vertexList.get(j);
        he.head.he = he; // make sure the vertex has at least one half edge that points to it.
        he.twin = halfEdges.get( twin );
        if ( he.twin != null ) he.twin.twin = he;
        halfEdges.put( p, he );        
        return he;        
    }    
    
    /** 
     * Reset the solutions for heat and distance
     */
    public void resetHeatAndDistanceSolution() {
    	for ( Vertex v : vertices ) {
    		v.u0 = v.constrained? 1 : 0;
    		v.ut = v.u0;
    		v.phi = 0;
    	}
    }
    
    /** 
     * Perform a specified number of projected Gauss-Seidel steps of the heat diffusion equation.
     * The current ut values stored in the vertices will be refined by this call.
     * @param GSSteps number of steps to take
     * @param t solution time
     */
    public void solveHeatFlowStep( int GSSteps, double t ) {    	
    	// Solve (A - t LC) u_t = u_0 with constrained vertices holding their ut value fixed
    	// Note that this is a backward Euler step of the heat diffusion.
    	for ( int i = 0; i < GSSteps; i++ ) {
    		for ( Vertex v : vertices ) {
    			if ( v.constrained ) continue;  // do nothing for the constrained vertex!
    			
    			// TODO: 7 write inner loop code for the PGS heat solve
    			double M_ii = v.area - (v.LCii * t);
    			
    			//clockwise order adjacent vertices
    			Vertex[] adjVertices = new Vertex[v.valence()];
    			int index = 0;
    			HalfEdge start = v.he;
        		HalfEdge loop = v.he;
        		do {
        			adjVertices[index] = loop.twin.head;
        			index++;
        			loop = loop.next.twin;
        		}while(!loop.equals(start));
    			
    			
        		
    			double sum = 0;
    			for (int j=0; j<v.valence(); j++) {
    				double M_ij = -t * v.LCij[j];
    				sum = sum + M_ij * adjVertices[j].ut;
    			}
    			
    			
    			
    			v.ut = (v.u0 - sum) / M_ii;	
    			
    			
    		}	
    	}
    }
    
    
    
    
    /**
     * Compute the gradient of heat at each face
     */
    public void updateGradu() {
    	// TODO: 8 update the gradient of u from the heat values, i.e., f.gradu for each Face f
    	
    	for (Face f: faces) {
    		HalfEdge start = f.he;
    		HalfEdge loop = f.he;
    		Vector3d gradient = new Vector3d();
    		do {
    		Vertex oppositeVertex = loop.next.head;
    		Vector3d v = new Vector3d(loop.head.p.x - loop.twin.head.p.x,
    				loop.head.p.y - loop.twin.head.p.y,
    				loop.head.p.z - loop.twin.head.p.z);
    		Vector3d product = new Vector3d();
    		product.cross(v, f.n);
    		
    		product = new Vector3d (product.x * oppositeVertex.ut,
    				product.y * oppositeVertex.ut,
    				product.z * oppositeVertex.ut);
    		
    		gradient.add(product);
    		
    		loop = loop.next;
    		}while (!loop.equals(start));
		
    		gradient.normalize();		
    		f.gradu = gradient;
    	}
    	
    	
    	
    }
    
    /** 
     * Compute the divergence of normalized gradients at the vertices
     */
    public void updateDivx() {
    	// TODO: 9 Update the divergence of the normalized gradients, ie., v.divX for each Vertex v
    	for (Vertex v: vertices) {
    		v.divX = 0;
    		double sum = 0;
    		
    		HalfEdge loop = v.he;
    		do {
    			Vector3d e1 = new Vector3d(loop.twin.head.p.x - v.p.x, 
    					loop.twin.head.p.y - v.p.y, 
    					loop.twin.head.p.z - v.p.z);
    			
    			Vector3d e2 = new Vector3d(loop.twin.next.head.p.x - v.p.x,
    					loop.twin.next.head.p.y - v.p.y,
    					loop.twin.next.head.p.z - v.p.z);
    			
    			Face f = loop.twin.leftFace;    			
    			double dot1 = e1.dot(f.gradu);
    			double dot2 = e2.dot(f.gradu);
    			
    			double cot_theta1 = 1.0 / Math.tan(angleWithNext (loop.twin.next));
        		double cot_theta2 = 1.0 / Math.tan(angleWithNext (loop.twin));
        		
        		sum += (cot_theta1*dot1) + (cot_theta2*dot2);
    			
    			loop = loop.twin.next.next;
    			
    		}while (!loop.equals(v.he));
    		
    		v.divX = sum / 2.0;
    		//System.out.println(v.divX);
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		/*
    		
    		v.divX = 0;
    		HalfEdge start = v.he.twin;
    		HalfEdge loop = v.he.twin;
    		double sum = 0;
    		do {
    		
    		Vector3d v1 = new Vector3d (loop.head.p.x - v.p.x,
    				loop.head.p.y - v.p.y,
    				loop.head.p.z - v.p.z);
    		Vector3d v2 =  new Vector3d (loop.next.next.twin.head.p.x - v.p.x,
    				loop.next.next.twin.head.p.y - v.p.y,
    				loop.next.next.twin.head.p.z - v.p.z
    				);
    		
    		Vector3d grad = loop.leftFace.gradu;
    		double dot1 = v1.dot(grad);
    		double dot2 = v2.dot(grad);
    		
    		double cot_theta1 = 1.0 / Math.tan(angleWithNext (loop.next));
    		double cot_theta2 = 1.0 / Math.tan(angleWithNext (loop));
    		
    		sum = sum + (cot_theta1 * dot1) + (cot_theta2 * dot2);
    		loop = loop.next.next.twin;
    		} while (!loop.equals(start));
    		
    		
    		v.divX = sum / 2.0;
    		System.out.println(v.divX);
    		*/
    		
    		
    	}
    	
    	
    	
    }
    
    /** Keep track of the maximum distance for debugging and colour map selection */
    double maxphi = 0 ;

    /**
     * Solves the distances
     * Uses Poisson equation, Laplacian of distance equal to divergence of normalized heat gradients.
     * This is step III in Algorithm 1 of the Geodesics in Heat paper, but here is done iteratively 
     * with a Gauss-Seidel solve of some number of steps to refine the solution whenever this method 
     * is called.
     * @param GSSteps number of Gauss-Seidel steps to take
     */
    public void solveDistanceStep( int GSSteps ) {		
    	for ( int i = 0; i < GSSteps; i++ ) {
    		for ( Vertex v : vertices ) {
    			// TODO: 10 Implement the inner loop of the Gauss-Seidel solve to compute the distances to each vertex, phi
    			
    			//LC phi = divx
    			if (Double.isNaN(v.phi)) 
    			{  
    				v.phi = 0;
    				continue;
    			}
    			   	
    			//all adjacent vertices   			
    			//clockwise order adjacent vertices
    			Vertex[] adjVertices = new Vertex[v.valence()];
    			int index = 0;
    			HalfEdge start = v.he;
        		HalfEdge loop = v.he;
        		do {
        			adjVertices[index] = loop.twin.head;
        			index++;
        			loop = loop.next.twin;
        		}while(loop != start);
        		
        		double sum = 0;
    			for (int j=0; j<v.valence(); j++) {
    				sum = sum + (v.LCij[j] * adjVertices[j].phi);
    			}
    			
    			v.phi = (v.divX - sum) / (v.LCii );
	
    		}    		
    	}
    	
    	// Note that the solution to step III is unique only up to an additive constant,
    	// final values simply need to be shifted such that the smallest distance is zero. 
    	// We also identify the max phi value here to identify the maximum geodesic and to 
    	// use adjusting the colour map for rendering
    	double minphi = Double.MAX_VALUE;
    	maxphi = Double.MIN_VALUE;
		for ( Vertex v : vertices ) {
			if ( v.phi < minphi ) minphi = v.phi;
			if ( v.phi > maxphi ) maxphi = v.phi;
		}	
		maxphi -= minphi;
		for ( Vertex v : vertices ) {
			v.phi -= minphi;
		}
    }
    
   
    /**
     * Computes the cotangent Laplacian weights at each vertex.
	 * You can assume no boundaries and a triangular mesh! 
	 * You should store these weights in an array at each vertex,
	 * and likewise store the associated "vertex area", i.e., 1/3 of
	 * the surrounding triangles and NOT scale your Laplacian weights
	 * by the vertex area (see heat solve objective requirements).
     */
    public void computeLaplacian() {
    	for ( Vertex v : vertices ) {
    		// TODO: 6 Compute the Laplacian and store as vertex weights, and cotan operator diagonal LCii and off diagonal LCij terms.
    		v.area = 0;
    		//the diagonal elements of the cotangent operator should be placed into LCii
    		v.LCii = 0; 
    		//the weights for the adjacent vertices placed in LCij
    		v.LCij = new double[ v.valence() ];
    		ArrayList<HalfEdge> adjEdges = new ArrayList<HalfEdge> ();
    		// find all adjacent edges ji
    		// clockwise order
    		HalfEdge start = v.he;
    		HalfEdge loop = v.he;
    		do {
    			adjEdges.add(loop);
    			loop = loop.next.twin;
    		}while(!loop.equals(start));
    		
    		double diagonal = 0; //ui
    		double offDiagonal = 0; //uj
    		int index = 0;
    		double vertexArea = 0;
    		//for each adjacent ji edge
    		for (HalfEdge adj : adjEdges) {
    			vertexArea = vertexArea + adj.leftFace.area;
    		}
    		v.area =  vertexArea/3.0;
    		
    		
    		
    		
    		for (HalfEdge adj : adjEdges) {
    			//diagonol terms
    			double alpha_ij = angleWithNext (adj.next);
    			double beta_ij = angleWithNext (adj.twin.next);
    			double cotan_alpha_ij = 1.0 / Math.tan(alpha_ij);
    			double cotan_beta_ij = 1.0 / Math.tan(beta_ij);
    			
    			diagonal = diagonal + (cotan_alpha_ij + cotan_beta_ij);
    			
    			offDiagonal = 0.5  * (cotan_alpha_ij + cotan_beta_ij);
    			v.LCij[index] = offDiagonal;
    			index++;
    			
    		}
    		diagonal = -0.5  * diagonal;
    		v.LCii = diagonal;
    		
    		
    		
    		
    		
    		
    	}
    }
    
    /** 
     * Computes the angle between the provided half edge and the next half edge
     * @param he specify which half edge
     * @return the angle in radians
     */
    private double angleWithNext( HalfEdge he ) {
    	// TODO: 6 Implement this function to compute the angle with next edge... you'll want to use this in a few places
    	
    	Vector3d ab = new Vector3d(he.head.p.x - he.twin.head.p.x,
    			he.head.p.y- he.twin.head.p.y,
    			he.head.p.z- he.twin.head.p.z);
    	ab.normalize();
    	
    	Vector3d cb = new Vector3d(he.head.p.x - he.next.head.p.x,
    			he.head.p.y- he.next.head.p.y,
    			he.head.p.z- he.next.head.p.z);
    	cb.normalize();
    	
    	double cosAlpha = ab.dot(cb);
    	double alpha = Math.acos(cosAlpha);
    	
    	return alpha;
    }
    
    public void shortestPath (Vertex v1, Vertex v2) {
    	
    }
    
    /**
     * Legacy drawing code for the half edge data structure by drawing each of its faces.
     * Legacy in that this code uses immediate mode OpenGL.  Per vertex normals are used
     * to draw the smooth surface if they are set in the vertices. 
     * @param drawable
     */
    public void display( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        for ( Face face : faces ) {
            HalfEdge he = face.he;
            if ( he.head.n == null ) { // don't have per vertex normals? use the face
                gl.glBegin( GL2.GL_POLYGON );
                Vector3d n = he.leftFace.n;
                gl.glNormal3d( n.x, n.y, n.z );
                HalfEdge e = he;
                do {
                	Point3d p = e.head.p;
                    gl.glVertex3d( p.x, p.y, p.z );
                    e = e.next;
                } while ( e != he );
                gl.glEnd();
            } else {
                gl.glBegin( GL2.GL_POLYGON );                
                HalfEdge e = he;
                do {
                	Point3d p = e.head.p;
                    Vector3d n = e.head.n;
                    gl.glNormal3d( n.x, n.y, n.z );
                    gl.glVertex3d( p.x, p.y, p.z );
                    e = e.next;
                } while ( e != he );
                gl.glEnd();
            }
        }
    }
    
}
