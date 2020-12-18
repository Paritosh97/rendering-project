#version 430
#define M_PI 3.14159265358979323846
#define MAX_SCENE_BOUNDS    10.0

struct DirLight
{   
    float ambient;
    float diffuse;
};
uniform DirLight dirLight;

uniform mat4 mat_inverse;
uniform mat4 persp_inverse;
uniform sampler2D envMap;
uniform vec3 center;
uniform float radius;
uniform bool transparent;
uniform float shininess;
uniform float eta;
uniform float eta_k;
uniform bool multipleSpheres;
uniform vec3 lightPosition;
uniform float lightIntensity;
uniform vec4 modelColor;
uniform int modelChoice;

in vec4 position;
in vec4 vertColor;
in vec4 vertNormal;
in vec2 textCoords;

out vec4 fragColor;

struct intersection_t
{
    vec3 start;
    vec4 color;
    vec3 normal;
    vec3 position;
    vec3 incomingDir;
    int sphere_id;
    float lambda;
    float fresnel;
};

struct sphere_t
{
    vec3 center;
    float radius;
    vec4 color;
};

float CalculateFresnel(float cos_td)
{
    float f_s = (pow(eta, 2) + pow(eta_k, 2) - (2*eta_k*cos_td) + pow(cos_td, 2)) / (pow(eta, 2) + pow(eta_k, 2) + (2*eta_k*cos_td) + pow(cos_td, 2));
    float f_p = ((pow(eta, 2) + pow(eta_k, 2))*pow(cos_td, 2) - (2*eta_k*cos_td) + 1) / ((pow(eta, 2) + pow(eta_k, 2))*pow(cos_td, 2) + (2*eta_k*cos_td) + 1);
    
    return (f_s + f_p)/2;
}

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

float chiGGX(float v)
{
    return v>0 ? 1.0:0.0;
}

/*
 * Microfacet normal distribution
 */
float D(float alpha, float cos_theta)
{
    float alpha2 = alpha * alpha;
    float cos2 = cos_theta * cos_theta;
    return (chiGGX(cos_theta) * alpha2) / (M_PI * pow(cos2 * alpha2 + (1 - cos2), 2));
}

/*
 * Geometric shadowing
 */
float G1(float alpha, float cos_theta)
{
    float cos2 = cos_theta*cos_theta;
    return 2 / (1 + sqrt(1 + alpha*alpha + (1 - cos2)/cos2));
}

float DistributionGGX(vec3 N, vec3 H)
{
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;
	
    float nom    = shininess*shininess;
    float denom  = (NdotH2 * (shininess*shininess - 1.0) + 1.0);
    denom        = M_PI * denom * denom;
	
    return nom / denom;
}

float CalculateG(float cos_t, float sin_t)
{
    float tan_t = sin_t/cos_t;
    return 2/(1+sqrt(1+(shininess*shininess*tan_t*tan_t)));
}

float CalcCookTorranceSpecular(float F, vec3 vertNormal, vec3 lightVector,  vec3 halfwayDir, vec3 eyeVector)
{
    float cos_ti = dot(vertNormal, lightVector);
    float sin_ti = length(cross(vertNormal, lightVector));
    float cos_to = dot(vertNormal, eyeVector);
    float sin_to = length(cross(vertNormal, eyeVector));
    
    float D = DistributionGGX(vertNormal, halfwayDir);    
    float G1_ti = CalculateG(cos_ti, sin_ti);
    float G1_td = CalculateG(cos_to, sin_to);

    float num = F * D * G1_ti * G1_td;
    float denom = 4 * cos_ti * cos_to;

    return num/denom;
}

vec4 getColorFromEnvironment(in vec3 direction)
{
    float u,v;
    u = (atan(-direction.x, direction.z) + M_PI) / (2 * M_PI);
    v = (asin(direction.y) + M_PI / 2) / M_PI;
    return texture(envMap, vec2(u, v));
}

bool raySphereIntersect(in vec3 start, in vec3 direction, out intersection_t intersection, sphere_t sphere)
{
    vec3 CP = start - sphere.center;
    float CP2 = length(CP) * length(CP);
    float r2 = sphere.radius * sphere.radius;
    float udotCP = dot(direction, CP);
    float udotCP2 = udotCP * udotCP;

    float b = 2.0 * udotCP;
    float c = CP2 - r2;
    float delta = pow(b, 2.0) - 4*c;
    if (delta >= 0)
    {
        float sqrtDelta = sqrt(delta);
        float lambda1 = (-b + sqrtDelta) / 2.0;
        float lambda2 = (-b - sqrtDelta) / 2.0;
        float lambda = min(lambda1, lambda2);
        intersection.position = start + lambda * direction;
        intersection.normal = normalize(intersection.position - sphere.center);
        intersection.incomingDir = direction;
        intersection.color = sphere.color;
        intersection.lambda = lambda;
        intersection.start = start;
        return true;
    }
    return false;
}

vec3 insideSphereIntersect(vec3 start, vec3 direction, vec3 center, float radius)
{
    vec3 CP = start - center;
    float CP2 = length(CP) * length(CP);
    float r2 = radius * radius;
    float udotCP = dot(direction, CP);
    float udotCP2 = udotCP * udotCP;

    float b = 2.0 * udotCP;
    float c = CP2 - r2;
    float delta = pow(b, 2.0) - 4*c;

    float sqrtDelta = sqrt(delta);
    float lambda1 = (-b + sqrtDelta) / 2.0;
    float lambda2 = (-b - sqrtDelta) / 2.0;
    float lambda = max(lambda1, lambda2);
    return start + lambda * direction;

}

bool closestSphereInterection(in sphere_t spheres[3], in vec3 start, in vec3 direction,in int id,  out intersection_t intersection_output)
{
    intersection_t intersection;
    intersection_t curr_intersection;

    intersection.lambda = MAX_SCENE_BOUNDS*100;
    bool intersected = false;

    for (int i = 0; i < 3; i++)
    {
        if(id == i)
        {
            continue;
        }
        if(raySphereIntersect(start, direction, curr_intersection, spheres[i]) && curr_intersection.lambda < intersection.lambda && curr_intersection.lambda > 0)
        {
            intersection = curr_intersection;
            intersected = true;
            intersection.sphere_id = i;
        }
    }

    intersection_output = intersection;
    return intersected;
}

vec4 RayBounceInside(vec3 u, vec3 intersection, int nb_bounces)
{
    vec3 normal = normalize(intersection - center);
    vec3 reflected_ray = normalize(reflect(u, normal));
    vec3 halfVector = normalize(reflected_ray + (-u));
    vec3 refracted_ray = normalize(refract(u, normal, eta));//

    float f = CalculateFresnel(dot(-u, halfVector));
    vec4 resultColor = f * getColorFromEnvironment(reflected_ray);
    float next_coef = 1 - f;
    vec3 incident_ray = refracted_ray;
    if(f < 0.999999)
    {
        for (int i = 0; i < 3; i++)
        {
            vec3 next_intersection = insideSphereIntersect(intersection, incident_ray, center, radius);
            normal = normalize(center - next_intersection);
            reflected_ray = normalize(reflect(incident_ray, normal));
            refracted_ray = normalize(refract(incident_ray, normal, 1/eta));
            halfVector = normalize((-incident_ray) + reflected_ray);

            f = CalculateFresnel(dot(-incident_ray, halfVector));
            resultColor += next_coef * (1-f) * getColorFromEnvironment(refracted_ray);

            next_coef = next_coef * f;
            incident_ray = reflected_ray;
            intersection = next_intersection;
        }
    }

    return resultColor;
}

vec4 GetColor(vec3 origin, vec4 inputColor, vec3 vertNormal, vec3 lightVector, vec3 halfwayDir, vec3 eyeVector, sphere_t spheres[3])
{
    vec4 ambient = dirLight.ambient * lightIntensity * inputColor;
    intersection_t intersection;
    if (closestSphereInterection(spheres, origin, lightVector, -1, intersection))
    {
        return ambient;
    }
    else
    {
        float cos_ti = dot(vertNormal, lightVector);
        float cos_td = dot(halfwayDir, lightVector);
        float F = CalcMetallicFresnel(cos_td);
        vec4 diffuse = dirLight.diffuse * max(cos_ti, 0) * inputColor * lightIntensity;
        vec4 specular = CalcCookTorranceSpecular(F, vertNormal, lightVector, halfwayDir, eyeVector) * lightIntensity * inputColor;
        return ambient + diffuse + specular;
    }
}

vec4 RayBounceOutside(sphere_t spheres[3], vec3 start, vec3 direction)
{
    intersection_t intersection;
    vec3 eyePosition = start;
    intersection_t intersection_stack[3];

    vec3 reflected_ray;
    vec3 lightVector;
    vec3 halfVector;
    vec3 eyeVector;
    vec4 resultColor = vec4(0.0, 0.0, 0.0, 0.0);

    intersection.lambda = MAX_SCENE_BOUNDS*100;
    int bounces = 0;
    int id_last_intersection = 0;
    vec3 incidentDir = direction;
    vec3 incidentPos = start;
    int id_sphere = -1;

    for (int i = 0; i < 3; i++)
    {
        if (closestSphereInterection(spheres, incidentPos, incidentDir, id_sphere, intersection))
        {
            intersection_stack[i] = intersection;
            id_sphere = intersection.sphere_id;
        }
        else
        {
            if(i == 0)
            {
                return vec4(0,0,0,0);
            }
            intersection_stack[i].color = vec4(0,0,0,0);
            id_last_intersection = i;
            break;
        }
        reflected_ray = normalize(reflect(intersection_stack[i].incomingDir, intersection_stack[i].normal));
        halfVector = normalize((-intersection_stack[i].incomingDir) + reflected_ray);
        if(i != 0)
        {
            intersection_stack[i].fresnel = CalcMetallicFresnel(dot((-intersection_stack[i].incomingDir), halfVector)) * intersection_stack[i-1].fresnel;
        }
        else
        {
            intersection_stack[i].fresnel = CalcMetallicFresnel(dot((-intersection_stack[i].incomingDir), halfVector));
        }
        incidentPos = intersection_stack[i].position;
        incidentDir = reflected_ray;
        id_last_intersection = i;
    }

    for (int i = id_last_intersection; i >= 0; i--)
    {
        reflected_ray = normalize(reflect(intersection_stack[i].incomingDir, intersection_stack[i].normal));
        lightVector = normalize(lightPosition - intersection_stack[i].position);
        eyeVector = normalize(-intersection_stack[i].incomingDir);
        halfVector = normalize(lightVector + eyeVector);
        if (i == id_last_intersection)
        {
            intersection_stack[i].color = dirLight.ambient * lightIntensity * intersection_stack[i].color;
        }
        else
        {
            intersection_stack[i].color = GetColor(intersection_stack[i].start, intersection_stack[i].color, intersection_stack[i].normal, lightVector, halfVector, eyeVector, spheres) + intersection_stack[i+1].color * intersection_stack[i].fresnel ;
        }
    }
    return intersection_stack[0].color;
}

void main(void)
{
    // Step 1: I need pixel coordinates. Division by w?
    vec4 worldPos = position;
    worldPos.z = 1; // near clipping plane
    worldPos = persp_inverse * worldPos;
    worldPos /= worldPos.w;
    worldPos.w = 0;
    worldPos = normalize(worldPos);

    // Step 2: ray direction:
    vec3 u = normalize((mat_inverse * worldPos).xyz);
    vec3 eye = (mat_inverse * vec4(0, 0, 0, 1)).xyz;

    if(multipleSpheres)
    {
        sphere_t spheres[3];
        spheres[0].center = vec3(1, 1, 1);
        spheres[0].radius = 0.4;
        spheres[0].color = vec4(1, 0, 0, 1);
        spheres[1].center = vec3(1, 0, 1);
        spheres[1].radius = 0.4;
        spheres[1].color = vec4(0, 1, 0, 1);
        spheres[2].center = vec3(0, 1, 1);
        spheres[2].radius = 0.3;
        spheres[2].color = vec4(0, 0, 1, 1);

        fragColor = RayBounceOutside(spheres, eye, u);
    } 
    else
    {
        vec4 resultColor;
        intersection_t intersection;
        sphere_t sphere;
        sphere.center = center;
        sphere.radius = radius;
        if (raySphereIntersect(eye, u, intersection, sphere))
        {
            if (transparent)
            {
                resultColor = RayBounceInside(u, intersection.position, 3);
            }
            else
            {
                vec3 normal = normalize(intersection.position - center);
                vec3 reflected_ray = normalize(reflect(u, normal));
                vec3 halfVector = normalize(reflected_ray + (-u));
                resultColor = CalcMetallicFresnel(dot(-u, halfVector)) * getColorFromEnvironment(reflected_ray);
            }
        }
        else
        {
            resultColor = getColorFromEnvironment(u);
        }
        fragColor = resultColor;
    }
}
