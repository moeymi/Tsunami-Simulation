// This shader visuzlizes the normal vector values on the mesh.
Shader "Example/URPUnlitShaderNormal"
{
    Properties
    {
        _BaseColor("Base Color", Color) = (1, 1, 1, 1)
        _MainTex("Albedo (RGB)", 2D) = "white" {}
        _Glossiness("Smoothness", Range(0,1)) = 0.5
        _Metallic("Metallic", Range(0,1)) = 0.0
        _Size("Particle Size", float) = 1
    }

        SubShader
    {
        Tags { "RenderType" = "Opaque" "RenderPipeline" = "UniversalRenderPipeline" }

        Pass
        {

            HLSLPROGRAM
            #pragma vertex vert            
            #pragma fragment frag
            #pragma multi_compile_instancing

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"            

            struct Attributes
            {
                float4 position   : POSITION;
                half3 normal        : NORMAL;
                uint instanceID : SV_InstanceID;
            };

            float _Size;
            half4 _BaseColor;
            StructuredBuffer<float3> _particlesBuffer;
            struct Varyings
            {
                float4 position : SV_POSITION;
                half3 normal : TEXCOORD0;
            };

            Varyings vert(Attributes IN)
            {
                float4 data = float4(_particlesBuffer[IN.instanceID], 0);

                float3 localPosition = IN.position.xyz;
                localPosition *= _Size;

                float3 worldPosition = localPosition + data;
                float4 finalPos = mul(UNITY_MATRIX_VP, float4(worldPosition, 1.0f));

                float3 worldNormal = IN.normal;

                Varyings OUT;
                OUT.position = finalPos;
                OUT.normal = TransformObjectToWorldNormal(IN.normal);
                return OUT;
            }

            half4 frag(Varyings IN) : SV_Target
            {
                return _BaseColor;
            }
            ENDHLSL
        }
    }
}