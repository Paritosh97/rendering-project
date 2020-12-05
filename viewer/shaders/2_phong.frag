#version 410

struct DirLight
{   
    float ambient;
    float diffuse;
};
uniform DirLight dirLight;

uniform float lightIntensity;
uniform bool blinnPhong;
uniform float shininess;
uniform float eta;
uniform float eta_k;
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
    vec4 halfwayDir = normalize(eyeVector + lightVector);

    // diffuse shading
    float diff = max(dot(normal, lightVector), 0.0);

    // specular shading
    vec4 reflectDir = reflect(-lightVector, normal);
    float cos_td = dot(halfwayDir.xyz, lightVector.xyz);
    float f_s = (pow(eta, 2) + pow(eta_k, 2) - (2*eta_k*cos_td) + pow(cos_td, 2)) / (pow(eta, 2) + pow(eta_k, 2) + (2*eta_k*cos_td) + pow(cos_td, 2));
    float f_p = ((pow(eta, 2) + pow(eta_k, 2))*pow(cos_td, 2) - (2*eta_k*cos_td) + 1) / ((pow(eta, 2) + pow(eta_k, 2))*pow(cos_td, 2) + (2*eta_k*cos_td) + 1);
    float fresnelCoefficient = (f_s + f_p)/2;
    float spec = pow(max(dot(normal, halfwayDir), 0.0), shininess);

    // combine results
    vec4 ambient  = dirLight.ambient  * vertColor;
    vec4 diffuse  = dirLight.diffuse  * diff * vertColor;
    vec4 specular = fresnelCoefficient * spec * vertColor;

    return (ambient + diffuse + specular);
}

void main( void )
{
    fragColor = lightIntensity * CalcBlinnPhong();
}
