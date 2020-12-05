#version 410

uniform mat4 matrix;
uniform mat4 perspective;
uniform mat3 normalMatrix;
uniform bool noColor;
uniform vec3 lightPosition;

// World coordinates
in vec4 vertex;
in vec4 normal;
in vec4 color;
in vec2 texcoords;

// Camera-space coordinates
out vec4 eyeVector;
out vec4 lightVector;
out vec4 lightSpace;
out vec4 vertColor;
out vec4 vertNormal;
out vec2 textCoords;

void main( void )
{
    if (noColor) vertColor = vec4(0.2, 0.6, 0.7, 1.0 );
    else vertColor = color;
    vertNormal.xyz = normalize(normalMatrix * normal.xyz);
    vertNormal.w = 0.0;
    textCoords = texcoords;

    // TODO: compute eyeVector, lightVector. 
    lightVector.xyz = normalize(lightPosition - vertex.xyz);
    lightVector.w = 0;
    eyeVector = normalize(vec4(matrix[3]) - vertex);

    gl_Position = perspective * matrix * vertex;
}
