#version 130
uniform vec3 lightPosition;         //needs to be normalized
uniform vec3 sphereColor;
uniform float sphereRadius;

in vec2 sphereLocal;
in float interpolatedDepth;

out vec4 fragColor;

void main()
{
    float distance = length(sphereLocal);
    if (distance > 1.0) discard;

    vec3 normal = vec3(sphereLocal, sqrt(1.0 - distance * distance));
    float depth = interpolatedDepth - sphereRadius * normal.z;
    gl_FragDepth = depth;

    float lightingIntensity = clamp(dot(lightPosition, normal), 0.0, 1.0);
    vec3 finalSphereColor = 0.1 * sphereColor           // diffuse
            + 0.5 * sphereColor * lightingIntensity     // ambient
            + 0.5 * pow(lightingIntensity, 50.0);       // specular

    fragColor = vec4(finalSphereColor, 1.0);
//    fragColor = vec4(1.0);
}
