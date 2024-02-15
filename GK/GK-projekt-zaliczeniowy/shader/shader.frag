#version 420

#define MAX_NLIGHTS 8

in FVertex {
    vec3 Colour;
    vec3 Position;
    vec3 Normal;
} In;

out vec4 out_Colour;

uniform mat4 mm;
uniform mat4 vm;
uniform mat4 pm;
uniform mat4 vpm;
uniform vec4 eyepos;


uniform vec4 ls_position;
uniform vec3 ls_ambient;
uniform vec3 ls_direct;
uniform vec3 ls_attenuation;
uniform uint ls_mask = 0x1;
uniform uint nls = 1;


vec3 posDifference (vec4 p, vec3 pos, out float dist)
{
    vec3 v;

    if (p.w != 0.0) {
        v = p.xyz/p.w-pos.xyz;
        dist = sqrt (dot (v, v));
    }
    else
    v = p.xyz;
    return normalize (v);
} /*posDifference*/

float attFactor (vec3 att, float dist)
{
    return 1.0/(((att.z*dist)+att.y)*dist+att.x);
} /*attFactor*/

vec3 LambertLighting (void)
{
    vec3  normal, lv, vv, Colour;
    float d, dist;
    uint  i, mask;

    normal = normalize (In.Normal);
    vv = posDifference (eyepos, In.Position, dist);

    Colour = vec3(0.0);
    for (i = 0, mask = 0x00000001;  i < nls;  i++, mask <<= 1)
    if ((ls_mask & mask) != 0) {
        Colour += ls_ambient * In.Colour;
        lv = posDifference (ls_position, In.Position, dist);
        d = dot (lv, normal);
        if (dot (vv, normal) > 0.0) {
            if (d > 0.0) {
                if (ls_position.w != 0.0)
                d *= attFactor (ls_attenuation, dist);
                Colour += (d * ls_direct) * In.Colour;
            }
        }
        else {
            if (d < 0.0) {
                if (ls_position.w != 0.0)
                d *= attFactor (ls_attenuation, dist);
                Colour -= (d * ls_direct) * In.Colour;
            }
        }
    }
    return clamp (Colour, 0.0, 1.0);
} /*LambertLighting*/

#define AGamma(colour) pow (colour, vec3(256.0/563.0))

void main (void)
{
    out_Colour = vec4 (AGamma ( LambertLighting () ), 1.0);
} /*main*/
