#version 400 core

uniform mat4 MVP;

layout (location = 0) in vec2 pos;
layout (location = 1) in vec4 vColor;

out vec4 color;

void main()
{
    gl_Position = MVP * vec4(pos, 0.0, 1.0);
    color = vColor;
}
