#version 330 core

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;
in vec4 VertexColor;

uniform sampler2D texture0;
uniform bool useTexture;
uniform bool enableLighting;

// Lighting uniforms
uniform vec3 lightPos;
uniform vec3 lightColor;
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

        // Ambient - stronger base lighting
        vec3 ambient = ambientStrength * lightColor;

        // Diffuse with half-lambert for softer shadows
        // Half-lambert wraps lighting around the model more
        vec3 lightDir = normalize(lightPos - FragPos);
        float NdotL = dot(norm, lightDir);
        float halfLambert = NdotL * 0.5 + 0.5; // Remap from [-1,1] to [0,1]
        halfLambert = halfLambert * halfLambert; // Square for falloff
        vec3 diffuse = halfLambert * lightColor;

        // Simple rim lighting for better visibility
        vec3 viewDir = normalize(viewPos - FragPos);
        float rim = 1.0 - max(dot(viewDir, norm), 0.0);
        rim = smoothstep(0.6, 1.0, rim);
        vec3 rimLight = rim * 0.2 * lightColor;

        result = (ambient + diffuse + rimLight) * baseColor.rgb;
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
