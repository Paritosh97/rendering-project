#version 410

#define PI 3.1415

struct DirLight
{   
    float ambient;
    float diffuse;
};
uniform DirLight dirLight;

uniform float lightIntensity;
uniform int modelChoice;
uniform bool useTexture;
uniform float shininess;
uniform float eta;
uniform float eta_k;
uniform sampler2D shadowMap;
uniform vec4 modelColor;

in vec4 eyeVector;
in vec4 lightVector;
in vec4 vertColor;
in vec4 vertNormal;
in vec4 lightSpace;

out vec4 fragColor;

float n_min(float r)
{
    return (1-r)/(1+r);
}

float n_max(float r)
{
    return (1+sqrt(r))/(1-sqrt(r)); 
}

float get_n(float r,float g)
{
    return n_min(r)*g + (1-g)*n_max(r);
}

float get_k2(float r, float n)
{
    float nr = (n+1)*(n+1)*r-(n-1)*(n-1);
    return nr/(1-r);
}

float get_r(float n, float k)
{
    return ((n-1)*(n-1)+k*k)/((n+1)*(n+1)+k*k);
}

float get_g(float n, float k)
{
    float r = get_r(n,k);
    return (n_max(r)-n)/(n_max(r)-n_min(r));
}

float CalculateFresnel(float cos_td)
{
    float f_s = (pow(eta, 2) + pow(eta_k, 2) - (2*eta_k*cos_td) + pow(cos_td, 2)) / (pow(eta, 2) + pow(eta_k, 2) + (2*eta_k*cos_td) + pow(cos_td, 2));
    float f_p = ((pow(eta, 2) + pow(eta_k, 2))*pow(cos_td, 2) - (2*eta_k*cos_td) + 1) / ((pow(eta, 2) + pow(eta_k, 2))*pow(cos_td, 2) + (2*eta_k*cos_td) + 1);
    
    return (f_s + f_p)/2;
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

float CalcBlinnPhongSpecular(float F,vec4 halfwayDir)
{    
    float spec = pow(max(dot(vertNormal, halfwayDir), 0.0), shininess);

    return (F * spec);
}

float CalcCookTorranceSpecular(float F, vec4 halfwayDir)
{
    float cos_ti = dot(vertNormal, lightVector);
    float sin_ti = length(cross(vertNormal.xyz, lightVector.xyz));
    float cos_to = dot(vertNormal, eyeVector);
    float sin_to = length(cross(vertNormal.xyz, eyeVector.xyz));
    
    float D = DistributionGGX(vertNormal.xyz, halfwayDir.xyz);    
    float G1_ti = CalculateG(cos_ti, sin_ti);
    float G1_td = CalculateG(cos_to, sin_to);

    float num = F * D * G1_ti * G1_td;
    float denom = 4 * cos_ti * cos_to;

    return num/denom;
}

float CalcMetallicFresnel(float theta)
{
    //clamp parameters
    float _r = clamp(eta,0.0,0.99);
    //compute n and k
    float n = get_n(_r,eta_k);
    float k2 = get_k2(_r,n);

    float c = cos(theta);
    float rs_num = n*n + k2 - 2*n*c + c*c;
    float rs_den = n*n + k2 + 2*n*c + c*c;
    float rs = rs_num/rs_den;
    
    float rp_num = (n*n + k2)*c*c - 2*n*c + 1;
    float rp_den = (n*n + k2)*c*c + 2*n*c + 1;
    float rp = rp_num/rp_den;

    return 0.5*(rs+rp);
}

void main( void )
{   
    vec4 ambient, diffuse, specular;

    vec4 halfwayDir = normalize(eyeVector + lightVector);
    float cos_td = dot(halfwayDir.xyz, lightVector.xyz);
    float F = CalculateFresnel(cos_td);

    float diff = max(dot(vertNormal, lightVector), 0.0);

    if(useTexture)
    {
        // ambient shading
        ambient  = dirLight.ambient  * vertColor;

        // diffuse shading
        diffuse  = dirLight.diffuse  * diff * vertColor;

        switch(modelChoice)
        {
            case 0:
                specular = lightIntensity * CalcBlinnPhongSpecular(F, halfwayDir) * vertColor;
                break;
            case 1:
                specular = lightIntensity * CalcCookTorranceSpecular(F, halfwayDir) * vertColor;
                break;
            case 2:
                F = CalcMetallicFresnel(dot(vertNormal, eyeVector));
                specular = lightIntensity * CalcCookTorranceSpecular(F, halfwayDir) * vertColor;
                break;
        }
    }
    else
    {
        // ambient shading
        ambient  = dirLight.ambient  * modelColor;

        // diffuse shading
        diffuse  = dirLight.diffuse  * diff * modelColor;

        switch(modelChoice)
        {
            case 0:
                specular = lightIntensity * CalcBlinnPhongSpecular(F, halfwayDir) * modelColor;
                break;
            case 1:
                specular = lightIntensity * CalcCookTorranceSpecular(F, halfwayDir) * modelColor;
                break;
            case 2:
                F = CalcMetallicFresnel(dot(vertNormal, eyeVector));
                specular = lightIntensity * CalcCookTorranceSpecular(F, halfwayDir) * modelColor;
                break;
        }
    }

    // combine results
    fragColor = ambient + diffuse + specular;
}
