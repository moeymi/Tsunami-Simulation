using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEditor;
using UnityEngine;
using Random = UnityEngine.Random;

struct Tri
{
    public Vector3 p0;
    public Vector3 p1;
    public Vector3 p2;
    public Vector3 normal;
}

public class GPU_Particle_Manager : MonoBehaviour
{
    [Header("Particle properties")]
    public float radius = 1f;  // particle radius
    public Mesh particleMesh;
    public Material material;
    public float mass = 4f;
    public float gasConstant = 8.314f;
    public float restDensity = 9f;
    public float viscosityCoefficient = 2.5f;
    public float[] g = { 0.0f, -9.81f, 0.0f };
    public float damping = -0.37f;
    // ReSharper disable InconsistentNaming
    private float radius2;
    private float radius3;
    private float radius4;
    private float radius5;
    private float mass2;

    [Header("Simulation properties")]
    public bool TsunamiMode;
    public bool VolcanoMode;
    public int numberOfParticles = 50000;
    public int dimensions = 100;
    public int maximumParticlesPerCell = 500;

    public float floorHeight;

    public Color baseColor = Color.blue;
    public Color surfaceColor = Color.white;

    [Range(0.1f, 4f)]
    public float colorModifier = 0.8f;


    [Header("Debug information")]
    public float averageFPS;

    private Vector3[] _particles;
    private Vector3[] _prevParticles;
    private Vector3[] _duplicateParticles;

    private Tri[] _tris;
    private Color[] _colors;

    private int[] _neighbourList;
    private int[] _neighbourCollisionList;

    private uint[] _hashGrid;
    public uint[] _hashGridTracker;

    private uint[] _collisionHashGrid;
    public uint[] _collisionHashGridTracker;

    private float[] _densities;
    private float[] _pressures;
    private Vector3[] _velocities;
    private Vector3[] _forces;

    private Camera mainCamera;

    float k1;
    float k2;
    float k3;

    float timer; 

    private ComputeBuffer _particlesBuffer;
    private ComputeBuffer _prevParticlesBuffer;
    private ComputeBuffer _duplicateParticlesBuffer;

    private ComputeBuffer _trisBuffer;
    private ComputeBuffer _argsBuffer;
    private ComputeBuffer _neighbourListBuffer;
    private ComputeBuffer _neighbourTrackerBuffer;

    private ComputeBuffer _neighbourCollisionListBuffer;
    private ComputeBuffer _neighbourCollisionTrackerBuffer;

    private ComputeBuffer _hashGridBuffer;
    private ComputeBuffer _hashGridTrackerBuffer;

    private ComputeBuffer _collisionHashGridBuffer;
    private ComputeBuffer _collisionHashGridTrackerBuffer;

    private ComputeBuffer _densitiesBuffer;
    private ComputeBuffer _pressuresBuffer;
    private ComputeBuffer _velocitiesBuffer;
    private ComputeBuffer _forcesBuffer;
    private ComputeBuffer _colorsBuffer;
    private static readonly int ParticlesBufferProperty = Shader.PropertyToID("_particlesBuffer");
    private static readonly int ColorsBufferProperty = Shader.PropertyToID("_colorsBuffer");

    public ComputeShader computeShader;
    private int clearHashGridKernel;
    private int recalculateHashGridKernel;
    private int recalculateCollisionHashGridKernel;
    private int buildNeighbourListKernel;
    private int buildCollisionNeighbourListKernel;
    private int computeForcesKernel;
    private int computeCollisionsKernel;
    private int integrateKernel;

    public Transform collidersParent;

    [StructLayout(LayoutKind.Sequential, Size = 28)]
    private struct Position
    {
        public float x;
        public float y;
        public float z;
    }

    private void Awake()
    {
        mainCamera = Camera.main;
        mainCamera.transform.position = new Vector3(dimensions / 2, dimensions/1.8f, -dimensions/4f);
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
        mass2 = mass * mass;

        k1 = (315 * mass) / (64 * Mathf.PI * Mathf.Pow(radius, 9));
        k2 = -(45 * mass) / (Mathf.PI * Mathf.Pow(radius, 6));
        k3 = (45 * viscosityCoefficient * mass) / (Mathf.PI * Mathf.Pow(radius, 6));

        InitTris();
        RespawnParticles();
        FindKernels();
        InitComputeShader();
        InitComputeBuffers();

        timer = 0; 
    }

    #region Initialisation

    private void RespawnParticles()
    {
        _particles = new Vector3[numberOfParticles];
        _prevParticles = new Vector3[numberOfParticles];
        _duplicateParticles = new Vector3[numberOfParticles * 3];
        _colors = new Color[numberOfParticles * 3];
        _densities = new float[numberOfParticles];
        _pressures = new float[numberOfParticles];
        _velocities = new Vector3[numberOfParticles];
        _forces = new Vector3[numberOfParticles];

        int counter = 0;
        float x_start_offset = 0 + radius;
        float y_start_offset = 0 + radius + floorHeight;
        float z_start_offset = 0 + radius;
        float x_end_offset = dimensions - radius;
        float y_end_offset = dimensions - radius;
        float z_end_offset = dimensions - radius;

        while (counter < numberOfParticles)
        {
            for (float y = y_start_offset; y < y_end_offset; y += radius/2)
                for (float x = x_start_offset; x < x_end_offset; x += radius)
                    for (float z = z_start_offset; z < z_end_offset; z += radius)
                    {
                        Vector3 startPos = new Vector3(x, y, z) + Vector3.one * Random.Range(-0.5f, 0.5f);
                        _particles[counter] = new Vector3
                        (
                            startPos.x,
                            startPos.y,
                            startPos.z
                        );
                        _densities[counter] = -1f;
                        _pressures[counter] = 0.0f;
                        _forces[counter] = Vector3.zero;
                        _velocities[counter] = new Vector3(-1, 0, 0);

                        if (++counter == numberOfParticles)
                        {
                            return;
                        }
                    
                    }
        }
    }

    private void FindKernels()
    {
        clearHashGridKernel = computeShader.FindKernel("ClearHashGrid");
        recalculateHashGridKernel = computeShader.FindKernel("RecalculateHashGrid");
        recalculateCollisionHashGridKernel = computeShader.FindKernel("RecalculateCollisionHashGrid");
        buildNeighbourListKernel = computeShader.FindKernel("BuildNeighbourList");
        buildCollisionNeighbourListKernel = computeShader.FindKernel("BuildCollisionNeighbourList");
        computeForcesKernel = computeShader.FindKernel("ComputeForces");
        computeCollisionsKernel = computeShader.FindKernel("ComputeCollisions");
        integrateKernel = computeShader.FindKernel("Integrate");
    }

    private void InitComputeShader()
    {
        computeShader.SetFloat("CellSize", radius * 2); // Setting cell-size h to particle diameter.
        computeShader.SetInt("Dimensions", dimensions);
        computeShader.SetInt("maximumParticlesPerCell", maximumParticlesPerCell);
        computeShader.SetFloat("radius", radius);
        computeShader.SetFloat("radius2", radius2);
        computeShader.SetFloat("radius3", radius3);
        computeShader.SetFloat("radius4", radius4);
        computeShader.SetFloat("radius5", radius5);
        computeShader.SetFloat("mass", mass);
        computeShader.SetFloat("mass2", mass2);
        computeShader.SetFloat("gasConstant", gasConstant);
        computeShader.SetFloat("restDensity", restDensity);
        computeShader.SetFloat("viscosityCoefficient", viscosityCoefficient);
        computeShader.SetFloat("damping", damping);
        computeShader.SetFloats("g", g);
        computeShader.SetFloats("epsilon", Mathf.Epsilon);
        computeShader.SetFloats("pi", Mathf.PI);
        computeShader.SetFloats("k1", k1);
        computeShader.SetFloats("k2", k2);
        computeShader.SetFloats("k3", k3);
        computeShader.SetVector("baseColor", baseColor);
        computeShader.SetVector("surfaceColor", surfaceColor);
        computeShader.SetFloat("colorModifier", colorModifier);
        computeShader.SetBool ("TsunamiMode", TsunamiMode);
        computeShader.SetBool("VolcanoMode", VolcanoMode);
    }

    // for calculations
    Vector3[] vertices;
    int[] triangles;
    Vector3[] normals;
    Matrix4x4 localToWorld;
    private void InitTris()
    {
        MeshFilter[] colliders = collidersParent.GetComponentsInChildren<MeshFilter>();
        List<Tri> trisList = new List<Tri>();
        for (int i = 0; i < colliders.Length; i++)
        {
            vertices = colliders[i].mesh.vertices;
            normals = colliders[i].mesh.normals;
            triangles = colliders[i].mesh.triangles;
            localToWorld = colliders[i].transform.localToWorldMatrix;
            for (int j = 0; j < triangles.Length; j += 3)
            {
                Vector3 p0 = localToWorld.MultiplyPoint3x4(vertices[triangles[j]]);
                Vector3 p1 = localToWorld.MultiplyPoint3x4(vertices[triangles[j + 1]]);
                Vector3 p2 = localToWorld.MultiplyPoint3x4(vertices[triangles[j + 2]]);
                Vector3 p3 = (p2 + p1 + p0)/3;
                Tri tri = new Tri
                {
                    normal = colliders[i].transform.rotation * normals[triangles[j]],
                    p0 = p0,
                    p1 = p1,
                    p2 = p2
                };
                trisList.Add(tri);
            }
        }
        _tris = trisList.ToArray();
    }
    void InitComputeBuffers()
    {
        uint[] args = {
            particleMesh.GetIndexCount(0),
            (uint) numberOfParticles * 3,
            particleMesh.GetIndexStart(0),
            particleMesh.GetBaseVertex(0),
            0
        };
        _argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _argsBuffer.SetData(args);

        _particlesBuffer = new ComputeBuffer(numberOfParticles, sizeof(float) * 3);
        _particlesBuffer.SetData(_particles);

        _prevParticlesBuffer = new ComputeBuffer(numberOfParticles, sizeof(float) * 3);
        _prevParticlesBuffer.SetData(_prevParticles);

        _duplicateParticlesBuffer = new ComputeBuffer(numberOfParticles * 3, sizeof(float) * 3);
        _duplicateParticlesBuffer.SetData(_duplicateParticles);

        _trisBuffer = new ComputeBuffer(_tris.Length, sizeof(float) * 3 * 4);
        _trisBuffer.SetData(_tris);

        _colorsBuffer = new ComputeBuffer(numberOfParticles * 3, sizeof(float) * 4);
        _colorsBuffer.SetData(_colors);

        _neighbourList = new int[numberOfParticles * maximumParticlesPerCell * 8];
        _neighbourCollisionList = new int[numberOfParticles * maximumParticlesPerCell * 8];

        _hashGrid = new uint[dimensions * dimensions * dimensions * maximumParticlesPerCell];
        _hashGridTracker = new uint[dimensions * dimensions * dimensions];

        _collisionHashGrid = new uint[dimensions * dimensions * dimensions * maximumParticlesPerCell];
        _collisionHashGridTracker = new uint[dimensions * dimensions * dimensions];

        _neighbourListBuffer = new ComputeBuffer(numberOfParticles * maximumParticlesPerCell * 8, sizeof(int));
        _neighbourListBuffer.SetData(_neighbourList);
        _neighbourTrackerBuffer = new ComputeBuffer(numberOfParticles, sizeof(int));

        _neighbourCollisionListBuffer = new ComputeBuffer(maximumParticlesPerCell * numberOfParticles * 8, sizeof(int));
        _neighbourCollisionListBuffer.SetData(_neighbourCollisionList);
        _neighbourCollisionTrackerBuffer = new ComputeBuffer(numberOfParticles, sizeof(int));

        _hashGridBuffer = new ComputeBuffer(dimensions * dimensions * dimensions * maximumParticlesPerCell, sizeof(uint));
        _hashGridBuffer.SetData(_hashGrid);

        _hashGridTrackerBuffer = new ComputeBuffer(dimensions * dimensions * dimensions, sizeof(uint));
        _hashGridTrackerBuffer.SetData(_hashGridTracker);

        _collisionHashGridBuffer = new ComputeBuffer(dimensions * dimensions * dimensions * maximumParticlesPerCell, sizeof(uint));
        _collisionHashGridBuffer.SetData(_collisionHashGrid);

        _collisionHashGridTrackerBuffer = new ComputeBuffer(dimensions * dimensions * dimensions, sizeof(uint));
        _collisionHashGridTrackerBuffer.SetData(_collisionHashGridTracker);

        _densitiesBuffer = new ComputeBuffer(numberOfParticles, sizeof(float));
        _densitiesBuffer.SetData(_densities);

        _pressuresBuffer = new ComputeBuffer(numberOfParticles, sizeof(float));
        _pressuresBuffer.SetData(_pressures);

        _velocitiesBuffer = new ComputeBuffer(numberOfParticles, sizeof(float) * 3);
        _velocitiesBuffer.SetData(_velocities);
        _forcesBuffer = new ComputeBuffer(numberOfParticles, sizeof(float) * 3);
        _forcesBuffer.SetData(_forces);

        computeShader.SetBuffer(clearHashGridKernel, "_hashGridTracker", _hashGridTrackerBuffer);

        computeShader.SetBuffer(recalculateHashGridKernel, "_particles", _particlesBuffer);
        computeShader.SetBuffer(recalculateHashGridKernel, "_hashGrid", _hashGridBuffer);
        computeShader.SetBuffer(recalculateHashGridKernel, "_hashGridTracker", _hashGridTrackerBuffer);

        computeShader.SetBuffer(recalculateCollisionHashGridKernel, "_tris", _trisBuffer);
        computeShader.SetBuffer(recalculateCollisionHashGridKernel, "_collisionHashGrid", _collisionHashGridBuffer);
        computeShader.SetBuffer(recalculateCollisionHashGridKernel, "_collisionHashGridTracker", _collisionHashGridTrackerBuffer);

        computeShader.SetBuffer(buildNeighbourListKernel, "_particles", _particlesBuffer);
        computeShader.SetBuffer(buildNeighbourListKernel, "_hashGrid", _hashGridBuffer);
        computeShader.SetBuffer(buildNeighbourListKernel, "_hashGridTracker", _hashGridTrackerBuffer);
        computeShader.SetBuffer(buildNeighbourListKernel, "_neighbourList", _neighbourListBuffer);
        computeShader.SetBuffer(buildNeighbourListKernel, "_neighbourTracker", _neighbourTrackerBuffer);

        computeShader.SetBuffer(buildCollisionNeighbourListKernel, "_particles", _particlesBuffer);
        computeShader.SetBuffer(buildCollisionNeighbourListKernel, "_prevParticles", _prevParticlesBuffer);
        computeShader.SetBuffer(buildCollisionNeighbourListKernel, "_tris", _trisBuffer);
        computeShader.SetBuffer(buildCollisionNeighbourListKernel, "_collisionHashGrid", _collisionHashGridBuffer);
        computeShader.SetBuffer(buildCollisionNeighbourListKernel, "_collisionHashGridTracker", _collisionHashGridTrackerBuffer);
        computeShader.SetBuffer(buildCollisionNeighbourListKernel, "_neighbourCollisionList", _neighbourCollisionListBuffer);
        computeShader.SetBuffer(buildCollisionNeighbourListKernel, "_neighbourCollisionTracker", _neighbourCollisionTrackerBuffer);

        computeShader.SetBuffer(computeForcesKernel, "_neighbourTracker", _neighbourTrackerBuffer);
        computeShader.SetBuffer(computeForcesKernel, "_neighbourList", _neighbourListBuffer);
        computeShader.SetBuffer(computeForcesKernel, "_particles", _particlesBuffer);
        computeShader.SetBuffer(computeForcesKernel, "_densities", _densitiesBuffer);
        computeShader.SetBuffer(computeForcesKernel, "_pressures", _pressuresBuffer);
        computeShader.SetBuffer(computeForcesKernel, "_velocities", _velocitiesBuffer);
        computeShader.SetBuffer(computeForcesKernel, "_forces", _forcesBuffer);

        computeShader.SetBuffer(computeCollisionsKernel, "_forces", _forcesBuffer);
        computeShader.SetBuffer(computeCollisionsKernel, "_tris", _trisBuffer);
        computeShader.SetBuffer(computeCollisionsKernel, "_densities", _densitiesBuffer);
        computeShader.SetBuffer(computeCollisionsKernel, "_particles", _particlesBuffer);
        computeShader.SetBuffer(computeCollisionsKernel, "_velocities", _velocitiesBuffer);
        computeShader.SetBuffer(computeCollisionsKernel, "_neighbourCollisionList", _neighbourCollisionListBuffer);
        computeShader.SetBuffer(computeCollisionsKernel, "_neighbourCollisionTracker", _neighbourCollisionTrackerBuffer);

        computeShader.SetBuffer(integrateKernel, "_particles", _particlesBuffer);
        computeShader.SetBuffer(integrateKernel, "_prevParticles", _prevParticlesBuffer);
        computeShader.SetBuffer(integrateKernel, "_duplicateParticles", _duplicateParticlesBuffer);
        computeShader.SetBuffer(integrateKernel, "_densities", _densitiesBuffer);
        computeShader.SetBuffer(integrateKernel, "_colors", _colorsBuffer);
        computeShader.SetBuffer(integrateKernel, "_forces", _forcesBuffer);
        computeShader.SetBuffer(integrateKernel, "_velocities", _velocitiesBuffer);
        computeShader.SetBuffer(integrateKernel, "_densities", _densitiesBuffer);
    }

    #endregion

    private void Start()
    {
        computeShader.Dispatch(recalculateCollisionHashGridKernel, _tris.Length, 1, 1);
    }
    Plane plane = new Plane(Vector3.up, 0);
    Vector3 worldPosition;
    void Update()
    {
        if (Input.GetMouseButtonDown(0))
        {
            float distance;
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if (plane.Raycast(ray, out distance))
            {
                worldPosition = ray.GetPoint(distance);
                VolcanoMode = true;
                timer = 0;
                computeShader.SetFloats("volcanoOrigin", new float[] { worldPosition.x, 0, worldPosition.z });
                computeShader.SetBool("VolcanoMode", VolcanoMode);
            }
        }
        if (VolcanoMode)
        {
            timer += Mathf.Min(0.08f, Time.deltaTime);
            if (timer > 0.2f)
            {
                timer = 0;
                VolcanoMode = false;
                computeShader.SetBool("VolcanoMode", VolcanoMode);
            }
        }

        if (TsunamiMode)
        {
            timer += Mathf.Min(0.08f, Time.deltaTime); 
            if (timer > 0.5)
            {
                timer = 0;
                TsunamiMode = false;
                computeShader.SetBool("TsunamiMode", TsunamiMode);
            }
        }

        computeShader.SetFloat("dt", Mathf.Min(0.08f, Time.deltaTime));
        computeShader.Dispatch(clearHashGridKernel, dimensions * dimensions * dimensions / 100, 1, 1);
        computeShader.Dispatch(recalculateHashGridKernel, numberOfParticles / 100, 1, 1);
        computeShader.Dispatch(buildNeighbourListKernel, numberOfParticles / 100, 1, 1);
        computeShader.Dispatch(buildCollisionNeighbourListKernel, numberOfParticles / 100, 1, 1);
        computeShader.Dispatch(computeForcesKernel, numberOfParticles / 100, 1, 1);
        computeShader.Dispatch(computeCollisionsKernel, numberOfParticles / 100, 1, 1);
        computeShader.Dispatch(integrateKernel, numberOfParticles / 100, 1, 1);
        material.SetBuffer(ParticlesBufferProperty, _duplicateParticlesBuffer);
        material.SetBuffer(ColorsBufferProperty, _colorsBuffer);
        Graphics.DrawMeshInstancedIndirect(particleMesh, 0, material, new Bounds(Vector3.zero, new Vector3(100.0f, 100.0f, 100.0f)), _argsBuffer, castShadows: 0);

        averageFPS = 1 / Time.deltaTime;
    }

    private void OnDestroy()
    {
        ReleaseBuffers();
    }

    private void ReleaseBuffers()
    {
        _particlesBuffer.Dispose();
        _prevParticlesBuffer.Dispose();
        _duplicateParticlesBuffer.Dispose();
        _trisBuffer.Dispose();
        _colorsBuffer.Dispose();
        _argsBuffer.Dispose();
        _neighbourListBuffer.Dispose();
        _neighbourTrackerBuffer.Dispose();
        _neighbourCollisionListBuffer.Dispose();
        _neighbourCollisionTrackerBuffer.Dispose();
        _hashGridBuffer.Dispose();
        _hashGridTrackerBuffer.Dispose();
        _collisionHashGridBuffer.Dispose();
        _collisionHashGridTrackerBuffer.Dispose();
        _densitiesBuffer.Dispose();
        _pressuresBuffer.Dispose();
        _velocitiesBuffer.Dispose();
        _forcesBuffer.Dispose();
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = Color.red;
        Gizmos.DrawWireCube(new Vector3(dimensions / 2f, dimensions / 2f, dimensions / 2f), Vector3.one * dimensions);
    }

    private void OnValidate()
    {
        computeShader?.SetVector("baseColor", baseColor);
        computeShader?.SetVector("surfaceColor", surfaceColor);
        computeShader?.SetFloat("colorModifier", colorModifier);
        computeShader?.SetBool("TsunamiMode", TsunamiMode);
        computeShader?.SetBool("VolcanoMode", VolcanoMode);
    }
}
