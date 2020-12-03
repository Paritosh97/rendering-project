#version 410

struct DirLight
{   
    float ambient;
    float diffuse;
    float specular;
};
uniform DirLight dirLight;

uniform float lightIntensity;
uniform bool blinnPhong;
uniform float shininess;
uniform float eta;
uniform sampler2D shadowMap;

in vec4 eyeVector;
in vec4 lightVector;
in vec4 vertColor;
in vec4 vertNormal;
in vec4 lightSpace;

out vec4 fragColor;

vec4 CalcBlinnPhong()
{
    vec4 normal = normalize(vertNormal);
    vec4 halfwayDir = normalize(lightVector);

    // diffuse shading
    float diff = max(dot(normal, lightVector), 0.0);

    // specular shading
    vec4 reflectDir = reflect(-lightVector, normal);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), shininess);

    // combine results
    vec4 ambient  = dirLight.ambient  * vertColor;
    vec4 diffuse  = dirLight.diffuse  * diff * vertColor;
    vec4 specular = dirLight.specular * spec * vertColor;

    return (ambient + diffuse + specular);
}

void main( void )
{
    fragColor = lightIntensity * CalcBlinnPhong();
}
