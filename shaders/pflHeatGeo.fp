//260761484
//Doreen He



#version 330 core

in vec3 camSpacePosition;
in vec3 camSpaceNormal;
in float utv;
in float phiv;

uniform vec3 lightCamSpacePosition;
uniform vec3 lightColor;
uniform vec3 materialDiffuse;

uniform int specularShininess;


float remap( float minval, float maxval, float curval ){
    return ( curval - minval ) / ( maxval - minval );
} 

void main(void) {
	
	vec3 v = normalize(-camSpacePosition);
	vec3 n = normalize(camSpaceNormal);
	vec3 l = normalize(lightCamSpacePosition - camSpacePosition);

	// TODO: 4, 11 Implement your GLSL per fragement lighting, heat colouring, and distance stripes here!
	
	

	const vec4 BLUE = vec4( 0.0, 0.0, 1.0, 1.0 );
	const vec4 RED   = vec4( 1.0, 0.0, 0.0, 1.0 );
	const vec4 BLACK   = vec4( 0.0, 0.0, 0.0, 1.0 );
	
	
	
	float u = clamp (utv, 0.0, 1.0);
	
	vec3 color ;
	if (u < 0.5)	
		color = (0.5 - utv) * 2 * vec3(0.0, 0.0, 1.0) + utv * 2 * vec3(0.6, 0.0, 0.0);
	else
		color = (utv - 0.5) * 2 * vec3(1.0, 1.0, 0.0) + (1 - utv) * 2 * vec3(0.6, 0.0, 0.0);
		
		
		
	//remap( 0.5, 0.5, utv )
		
	

	float diffuse = max(0.0, dot(n, l));
	vec3 halfAngle = normalize(v + l);
	float specular = max(0.0, dot(n, halfAngle));
	if (diffuse == 0.0)
		specular = 0;
	else
    	specular = pow(specular, specularShininess);
    	
    // lambertian and sprcularLight
	vec3 lambertian = materialDiffuse * diffuse;
    vec3 specularLight =  specular * lightColor;  
  	vec3 allLight =  lambertian + specularLight;
  	  
    
  
    float stripe = smoothstep(0.8, 0.85, mod(phiv* 10, 1.0)) * (1 - smoothstep(0.9, 0.95, mod(phiv*10, 1.0)));
    
    if (stripe==1.0) 
    	if (u<0.5)
    		color = vec3(0.0, 0.0, 1.0);
    	else
    		color = vec3(1.0, 0.0, 0.0);
    	
    
    
    gl_FragColor = vec4( allLight.x + color.x,  allLight.y + color.y , allLight.z + color.z, 1) ;
    
    
    
	// can use this to initially visualize the normal	
    // gl_FragColor = vec4( n.xyz * 0.5 + vec3( 0.5, 0.5,0.5 ), 1 );
    
    
}
