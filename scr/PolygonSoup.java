

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

/**
 * Simple implementation of a loader for a polygon soup
 */
public class PolygonSoup {

    /** List of vertex objects used in the mesh. */
    public ArrayList<Vertex> vertexList = new ArrayList<Vertex>();
    
    /** List of faces, where each face is a list of integer indices into the vertex list. */
    public ArrayList<int[]> faceList = new ArrayList<int[]>();
    
    /** Map for keeping track of how many n-gons we have for each n */
    private TreeMap<Integer,Integer> faceSidesHistogram = new TreeMap<Integer,Integer>();
    
    /** A string summarizing the face sides histogram */
    public String soupStatistics;
    
    /**
     * Creates a polygon soup by loading an OBJ file
     * @param file
     */
    public PolygonSoup(String file) {
        try {
            FileInputStream fis = new FileInputStream( file );
            InputStreamReader isr = new InputStreamReader( fis );
            BufferedReader reader = new BufferedReader( isr );
            String line;
            while ((line = reader.readLine()) != null) {
                if ( line.startsWith("v ") ) {
                    parseVertex(line);
                } else if ( line.startsWith("f ") ) {
                    parseFace(line);
                } 
            }
            reader.close();
            isr.close();
            fis.close();
                        
            soupStatistics = file + "\n" + "faces = " +faceList.size() + "\nverts = " + vertexList.size() + "\n";
            for ( Map.Entry<Integer,Integer> e : faceSidesHistogram.entrySet() ) {
                soupStatistics += e.getValue() + " ";
                if ( e.getKey() == 3 ) {
                    soupStatistics += "triangles\n";
                } else if ( e.getKey() == 4 ) {
                    soupStatistics += "quadrilaterals\n";
                } else {
                    soupStatistics += e.getKey() + "-gons\n";
                }
            }
            System.out.println( soupStatistics );
            
            //  compute a bounding box and scale and center the geometry
            /*
            Compute an axis aligned bounding box for the loaded mesh vertices, 
            and then modify the vertex positions so that they bounding box is centered at the origin 
            and has a length of 10 in the largest dimension. 
            Note that one could instead use glScale and glTranslate 
            to apply a modeling transform to the geometry, but applying the scale directly 
            to the mesh can help you to validate your geodesic calculations later 
            (e.g., for the sphere). Notice how small meshes, such as headtri.obj, 
            fill the screen nicely once you've completed this objective. 
            */
            
            // find bounding x,y,z limits
            double minX = vertexList.get(0).p.x;
            double maxX = minX;
            double minY = vertexList.get(0).p.y;
            double maxY = minY;
            double minZ = vertexList.get(0).p.z;
            double maxZ = minZ;
                    
            for (Vertex v : vertexList) {
            	if (v.p.x <= minX) minX = v.p.x;
            	if (v.p.x >= maxX) maxX = v.p.x;
            	if (v.p.y <= minY) minY = v.p.y;
            	if (v.p.y >= maxY) maxY = v.p.y;        
            	if (v.p.z <= minZ) minZ = v.p.z;            	
            	if (v.p.z >= maxZ) maxZ = v.p.z;
            }
            
            //center of the meshes
            Point3d center = new Point3d ( (minX + maxX)/2, (minY + maxY)/2, (minZ + maxZ)/2 );
            
            //scale size
            double lengthX = maxX - minX;
            double lengthY = maxY - minY;
            double lengthZ = maxZ - minZ;
            
            double maxLength = 0;
            if (lengthX > maxLength) maxLength = lengthX;
            if (lengthY > maxLength) maxLength = lengthY;
            if (lengthZ > maxLength) maxLength = lengthZ;
            
            double scale = 10/maxLength;
            
            //center and scale the meshes
            for (Vertex v: vertexList ) {            	
            	v.p.x = (v.p.x - center.x) * scale;
            	v.p.y = (v.p.y - center.y) * scale;
            	v.p.z = (v.p.z - center.z) * scale;
            }
            
                    
    
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Parses a vertex definition from a line in an obj file, and 
     * directly inserts it into the vertex list.
     * Assumes that there are three components.
     * @param newline
     * @return a new vertex object
     */
    private Vertex parseVertex(String newline) {        
        // Remove the tag "v "
        newline = newline.substring(2, newline.length());
        StringTokenizer st = new StringTokenizer(newline, " ");
        Vertex v = new Vertex();
        v.p.x = Double.parseDouble(st.nextToken());
        v.p.y = Double.parseDouble(st.nextToken());
        v.p.z = Double.parseDouble(st.nextToken());
        v.index = vertexList.size();
        vertexList.add( v );
        return v;
    }

    /**
     * Gets the list of indices for a face from a string in an obj file.
     * Simply ignores texture and normal information for simplicity
     * @param newline
     * @return list of indices
     */
    private int[] parseFace(String newline) {
        // Remove the tag "f "
        newline = newline.substring(2, newline.length());
        // vertex/texture/normal tuples are separated by spaces.
        StringTokenizer st = new StringTokenizer(newline, " ");
        int count = st.countTokens();
        int v[] = new int[count];
        for (int i = 0; i < count; i++) {
            // first token is vertex index... we'll ignore the rest (if it exists)
            StringTokenizer st2 = new StringTokenizer(st.nextToken(),"/");
            v[i] = Integer.parseInt(st2.nextToken()) - 1; // want zero indexed vertices!            
        }
        Integer n = faceSidesHistogram.get( count );
        if ( n == null ) {
            faceSidesHistogram.put( count, 1 );
        } else {
            faceSidesHistogram.put( count, n + 1 );
        }
        faceList.add( v );
        return v;
    }    

    /**
     * Draw the polygon soup using legacy immediate mode OpenGL
     * @param drawable
     */
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();
        // assume triangular faces!
        Vector3d v1 = new Vector3d();
        Vector3d v2 = new Vector3d();
        Vector3d n = new Vector3d();
        for ( int[] faceVertex : faceList ) {
            Point3d p0 = vertexList.get( faceVertex[0] ).p;
            Point3d p1 = vertexList.get( faceVertex[1] ).p;
            Point3d p2 = vertexList.get( faceVertex[2] ).p;
            v1.sub( p1,p0 );
            v2.sub( p2,p1 );
            n.cross( v1, v2 );
            gl.glBegin( GL2.GL_POLYGON );
            gl.glNormal3d( n.x, n.y, n.z );
            for ( int i = 0; i < faceVertex.length; i++ ) {
                Point3d p = vertexList.get( faceVertex[i] ).p;
                gl.glVertex3d( p.x, p.y, p.z );
            }
            gl.glEnd();
        }        
    }
}
