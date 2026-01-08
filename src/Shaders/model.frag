#version 330 core

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;
in vec4 VertexColor;

uniform sampler2D texture0;
uniform bool useTexture;
uniform bool enableLighting;

// Multi-light support (4 lights: Right, Left, Front, Rear)
#define MAX_LIGHTS 4
uniform vec3 lightPos[MAX_LIGHTS];
uniform vec3 lightColor[MAX_LIGHTS];
uniform bool lightEnabled[MAX_LIGHTS];

// View position for specular/rim
uniform vec3 viewPos;
uniform float ambientStrength;

// Alpha test
uniform bool enableAlphaTest;
uniform float alphaRef;

// Base alpha multiplier (set based on blend mode)
uniform float baseAlpha;

out vec4 FragColor;

void main()
{
    // Sample texture or use white if no texture
    vec4 texColor = useTexture ? texture(texture0, TexCoord) : vec4(1.0);

    // Combine with vertex color
    vec4 baseColor = texColor * VertexColor;

    // Alpha test
    if (enableAlphaTest && baseColor.a <= alphaRef)
        discard;

    vec3 result;

    if (enableLighting)
    {
        vec3 norm = normalize(Normal);
        vec3 viewDir = normalize(viewPos - FragPos);

        // Ambient base
        vec3 ambient = ambientStrength * vec3(1.0);
        vec3 totalDiffuse = vec3(0.0);

        // Accumulate lighting from all enabled lights
        for (int i = 0; i < MAX_LIGHTS; i++)
        {
            if (lightEnabled[i])
            {
                vec3 lightDir = normalize(lightPos[i] - FragPos);
                float NdotL = dot(norm, lightDir);

                // Half-lambert for softer shadows
                float halfLambert = NdotL * 0.5 + 0.5;
                halfLambert = halfLambert * halfLambert;

                totalDiffuse += halfLambert * lightColor[i];
            }
        }

        // Simple rim lighting for better visibility
        float rim = 1.0 - max(dot(viewDir, norm), 0.0);
        rim = smoothstep(0.6, 1.0, rim);
        vec3 rimLight = rim * 0.15 * vec3(1.0);

        result = (ambient + totalDiffuse + rimLight) * baseColor.rgb;
    }
    else
    {
        // When lighting is disabled, boost brightness slightly
        // to compensate for dark vertex colors common in FF7 models
        result = baseColor.rgb * 1.2;
    }

    // Clamp to valid range
    result = clamp(result, 0.0, 1.0);

    // Apply base alpha (set based on blend mode - 0.5 for BLEND_AVG, 1.0 otherwise)
    FragColor = vec4(result, baseColor.a * baseAlpha);
}
