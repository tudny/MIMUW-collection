#version 420

layout(location=0) in vec4 in_Position;
layout(location=1) in vec3 in_Colour;

out Vertex {
    vec3 Colour;
    vec3 Position;
} Out;


uniform mat4 mm;
uniform mat4 vm;
uniform mat4 pm;
uniform mat4 vpm;
uniform vec4 eyepos;


void main (void)
{
    vec4 Pos;

    Pos = mm * in_Position;
    Out.Position = Pos.xyz / Pos.w;
    Out.Colour = in_Colour;
    gl_Position = vpm * Pos;
} /*main*/
