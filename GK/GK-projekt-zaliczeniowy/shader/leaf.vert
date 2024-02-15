#version 420

#define MAX_INSTANCES 255

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 color;

uniform mat4 vpm;

out Vertex {
    vec3 Colour;
    vec3 Position;
} Out;

uniform mat4 modelMatrix[MAX_INSTANCES];

void main(void) {
    mat4 model = modelMatrix[gl_InstanceID];
    vec4 pos = model * vec4(position, 1.0);

    Out.Position = pos.xyz / pos.w;
    Out.Colour = color;
    gl_Position = vpm * pos;
}
