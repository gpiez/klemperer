#version 130
//position is always the same for all 4 vertices
//TODO use geometry shader
in vec4 position;

out vec2 sphereLocal;
out float interpolatedDepth;

uniform mat4 modelView;
uniform mat4 projection;
uniform float sphereRadius;

// DisplayWidget::initializeGL
const vec2 sphereLocalCoord[4] = vec2[4]( vec2(-1, -1), vec2(1, -1), vec2(-1, 1), vec2(1, 1) );

void main()
{
    vec4 displacedPos = modelView * position;
    sphereLocal = sphereLocalCoord[gl_VertexID];

    displacedPos.xy = displacedPos.xy + sphereLocal * vec2(sphereRadius);
    displacedPos = projection * displacedPos;

    gl_Position = displacedPos;
    interpolatedDepth = displacedPos.z;
}

