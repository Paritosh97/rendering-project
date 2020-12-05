#version 410

#define PI 3.1415

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

float CalculateFresnel(float cos_td)
{
    float f_s = (pow(eta, 2) + pow(eta_k, 2) - (2*eta_k*cos_td) + pow(cos_td, 2)) / (pow(eta, 2) + pow(eta_k, 2) + (2*eta_k*cos_td) + pow(cos_td, 2));
    float f_p = ((pow(eta, 2) + pow(eta_k, 2))*pow(cos_td, 2) - (2*eta_k*cos_td) + 1) / ((pow(eta, 2) + pow(eta_k, 2))*pow(cos_td, 2) + (2*eta_k*cos_td) + 1);
    
    return (f_s + f_p)/2;
}

vec4 CalcBlinnPhong()
{
    vec4 normal = normalize(vertNormal);
    vec4 halfwayDir = normalize(eyeVector + lightVector);

    // diffuse shading
    float diff = max(dot(normal, lightVector), 0.0);

    // specular shading
    vec4 reflectDir = reflect(-lightVector, normal);
    
    float cos_td = dot(halfwayDir.xyz, lightVector.xyz);
    float fresnelCoefficient = CalculateFresnel(cos_td);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), shininess);

    // combine results
    vec4 ambient  = dirLight.ambient  * vertColor;
    vec4 diffuse  = dirLight.diffuse  * diff * vertColor;
    vec4 specular = fresnelCoefficient * spec * vertColor;

    return (ambient + diffuse + specular);
}

float DistributionGGX(vec3 N, vec3 H)
{
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;
	
    float nom    = shininess*shininess;
    float denom  = (NdotH2 * (shininess*shininess - 1.0) + 1.0);
    denom        = PI * denom * denom;
	
    return nom / denom;
}

float CalculateG(float cos_t, float sin_t)
{
    float tan_t = sin_t/cos_t;
    return 2/(1+sqrt(1+(shininess*shininess*tan_t*tan_t)));
}

vec4 CalcCookTorrance()
{
    vec4 normal = normalize(vertNormal);
    vec4 halfwayDir = normalize(eyeVector + lightVector);

    // diffuse shading
    float diff = max(dot(normal, lightVector), 0.0);

    // specular shading
    vec4 reflectDir = reflect(-lightVector, normal);
    
    float cos_td = dot(halfwayDir.xyz, lightVector.xyz);
    float cos_ti = dot(normal, lightVector);
    float sin_ti = length(cross(normal.xyz, lightVector.xyz));
    float cos_to = dot(normal, eyeVector);
    float sin_to = length(cross(normal.xyz, eyeVector.xyz));
    
    float D = DistributionGGX(normal.xyz, halfwayDir.xyz);    
    float G1_ti = CalculateG(cos_ti, sin_ti);
    float G1_td = CalculateG(cos_to, sin_to);
    float F = CalculateFresnel(cos_td);

    float num = F * D * G1_ti * G1_td;
    float denom = 4 * cos_ti * cos_to;

    // combine results
    vec4 ambient  = dirLight.ambient  * vertColor;
    vec4 diffuse  = dirLight.diffuse  * diff * vertColor;
    vec4 specular = vertColor * num/denom;

    return (ambient + diffuse + specular);
}

void main( void )
{
    if(blinnPhong)
        fragColor = lightIntensity * CalcBlinnPhong();
    else
        fragColor = lightIntensity * CalcCookTorrance();
}
