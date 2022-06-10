using System.Collections.Generic;
using UnityEngine;
using Debug = UnityEngine.Debug;

public class SPH_Manager : MonoBehaviour
{
    #region Simulation Constants
    [Header("Particle properties")]
    public float radius = 1f; // particle radius, interaction radius h
    public GameObject particlePrefab;
    public float particleRadius = 1f;
    public float mass = 4f;
    public float viscosityCoefficient = 2.5f;
    public Vector3 g = new Vector3(0.0f, -9.81f, 0.0f);
    public float gasConstant = 1f;

    [SerializeField]
    private float restDensity = 1f;
    [SerializeField]
    private float damping = -0.5f;

    [Header("Simulation space properties")]
    public bool tsunamiMode = false;
    public int numberOfParticles = 1000;
    public int dimensions = 10;
    public int maximumParticlesPerCell = 500;

    [Header("Debug information")]
    [Tooltip("Tracks how many neighbours each particleIndex has in " + nameof(_neighbourList))]
    [SerializeField]
    private int[] _neighbourTracker;
    private int[] _neighbourCollisionTracker;

    public float averageFPS;
    #endregion

    private GameObject[] _particles;
    private int[] _neighbourList;
    private Triangle[] _neighbourCollisionList;
    private readonly Dictionary<int, List<int>> _hashGrid = new Dictionary<int, List<int>>();
    private readonly Dictionary<int, HashSet<Triangle>> _dynamicCollisionHashGrid = new Dictionary<int, HashSet<Triangle>>();
    private readonly Dictionary<int, HashSet<Triangle>> _staticCollisionHashGrid = new Dictionary<int, HashSet<Triangle>>();

    private float radius2;
    private float radius3;
    private float radius4;
    private float radius5;

    private float[] densities;
    private float[] pressures;
    private Vector3[] forces;
    private Vector3[] velocities;

    float k1;
    float k2;
    float k3;

    private MeshFilter[] staticColliders;
    private MeshFilter[] dynamicColliders;

    [SerializeField]
    private Transform staticCollidersParent;
    [SerializeField]
    private Transform dynamicCollidersParent;

    private void Awake()
    {
        RespawnParticles();
        InitNeighbourHashing();
        staticColliders = staticCollidersParent.GetComponentsInChildren<MeshFilter>();
        dynamicColliders = dynamicCollidersParent.GetComponentsInChildren<MeshFilter>();
        InitColliders(dynamicColliders, false);
        InitColliders(staticColliders, true);
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;

        k1 = (315 * mass) / (64 * Mathf.PI * Mathf.Pow(radius, 9));
        k2 = -(45 * mass) / (Mathf.PI * Mathf.Pow(radius, 6));
        k3 = (45 * viscosityCoefficient * mass) / (Mathf.PI * Mathf.Pow(radius, 6)); 
    }

    #region Initialisation

    private void InitColliders(MeshFilter[] colliders, bool isStatic)
    {

        // 2. Recalculate hashes of each mesh.
        for (int i = 0; i < colliders.Length; i++)
        {
            vertices = colliders[i].mesh.vertices;
            normals = colliders[i].mesh.normals;
            triangles = colliders[i].mesh.triangles;
            localToWorld = colliders[i].transform.localToWorldMatrix;
            for (int j = 0; j < triangles.Length; j += 3)
            {
                Vector3 facePos = (vertices[triangles[j]] + vertices[triangles[j + 1]] + vertices[triangles[j + 2]]) / 3;
                Vector3 p1 = localToWorld.MultiplyPoint3x4(vertices[triangles[j]]);
                Vector3 p2 = localToWorld.MultiplyPoint3x4(vertices[triangles[j + 1]]);
                Vector3 p3 = localToWorld.MultiplyPoint3x4(vertices[triangles[j + 2]]);
                float distance1 = Vector3.Distance(p1, p2);
                float distance2 = Vector3.Distance(p2, p3);
                float distance3 = Vector3.Distance(p1, p3);
                Triangle face = new Triangle(p1, p2, p3, colliders[i].transform.rotation * normals[triangles[j]]);
                if (distance1 > distance2 && distance1 > distance3)
                {
                    int t1 = Mathf.CeilToInt(distance2 / HashGrid.CellSize);
                    int t2 = Mathf.CeilToInt(distance3 / HashGrid.CellSize);
                    for (int h = 0; h <= t1; h++)
                    {
                        for (int l = 0; l <= t2; l++)
                        {
                            Vector3 potentialCell = p3 + ((p2 - p3).normalized * HashGrid.CellSize / 2 * h) + ((p1 - p3).normalized * HashGrid.CellSize / 2 * l);
                            if (face.PointInTriangle(potentialCell))
                            {
                                Debug.DrawLine(p1, potentialCell, Color.white, 5f);
                                if (!isStatic)
                                    _dynamicCollisionHashGrid[HashGrid.Hash(HashGrid.GetCell(potentialCell))].Add(face);
                                else
                                    _staticCollisionHashGrid[HashGrid.Hash(HashGrid.GetCell(potentialCell))].Add(face);
                            }
                        }
                    }
                }
                else if (distance2 > distance1 && distance2 > distance3)
                {
                    int t1 = Mathf.CeilToInt(distance1 / HashGrid.CellSize);
                    int t2 = Mathf.CeilToInt(distance3 / HashGrid.CellSize);
                    for (int h = 0; h <= t1; h++)
                    {
                        for (int l = 0; l <= t2; l++)
                        {
                            Vector3 potentialCell = p1 + ((p3 - p1).normalized * HashGrid.CellSize / 2 * l) + ((p2 - p1).normalized * HashGrid.CellSize / 2 * h);
                            if (face.PointInTriangle(potentialCell))
                            {
                                Debug.DrawLine(p1, potentialCell, Color.white, 5f);
                                if (!isStatic)
                                    _dynamicCollisionHashGrid[HashGrid.Hash(HashGrid.GetCell(potentialCell))].Add(face);
                                else
                                    _staticCollisionHashGrid[HashGrid.Hash(HashGrid.GetCell(potentialCell))].Add(face);
                            }
                        }
                    }
                }
                else
                {
                    int t1 = Mathf.CeilToInt(distance1 / HashGrid.CellSize);
                    int t2 = Mathf.CeilToInt(distance2 / HashGrid.CellSize);
                    int h = 0, l = 0;
                    for (h = 0; h <= t1; h++)
                    {
                        for (l = 0; l <= t2; l++)
                        {
                            Vector3 potentialCell = p2 + ((p3 - p2).normalized * HashGrid.CellSize / 2f * l) + ((p1 - p2).normalized * HashGrid.CellSize / 2 * h);
                            if (face.PointInTriangle(potentialCell))
                            {
                                Debug.DrawLine(p1, potentialCell, Color.white, 5f);
                                if (!isStatic)
                                    _dynamicCollisionHashGrid[HashGrid.Hash(HashGrid.GetCell(potentialCell))].Add(face);
                                else
                                    _staticCollisionHashGrid[HashGrid.Hash(HashGrid.GetCell(potentialCell))].Add(face);
                            }
                        }
                    }
                }
            }
        }
    }

    private void RespawnParticles()
    {
        _particles = new GameObject[numberOfParticles];
        densities = new float[numberOfParticles];
        pressures = new float[numberOfParticles];
        forces = new Vector3[numberOfParticles];
        velocities = new Vector3[numberOfParticles];

        int particlesPerDimension = Mathf.CeilToInt(Mathf.Pow(numberOfParticles, 1f / 3f));

        int counter = 0;

        GameObject parentParticle = new GameObject();
        float x_start_offset = 0 + particleRadius;
        float y_start_offset = 0 + particleRadius*5;
        float z_start_offset = 0 + particleRadius;
        float x_end_offset = dimensions - particleRadius;
        float y_end_offset = dimensions - particleRadius;
        float z_end_offset = dimensions - particleRadius;

        while (counter < numberOfParticles)
        {
            for (float y = y_start_offset; y < y_end_offset; y +=particleRadius)
                for (float x = x_start_offset; x < x_end_offset; x += particleRadius)
                    for (float z = z_start_offset; z < z_end_offset; z += particleRadius)
                    {
                        
                        Vector3 startPos = new Vector3(x,y,z) + Vector3.one * Random.Range(-0.5f, 0.5f);
                        _particles[counter] = Instantiate(particlePrefab);
                        _particles[counter].transform.parent = parentParticle.transform;
                        _particles[counter].transform.position = startPos;
                        _particles[counter].transform.localScale = new Vector3(particleRadius, particleRadius, particleRadius);

                        densities[counter] = -1f;
                        pressures[counter] = 0.0f;
                        forces[counter] = Vector3.zero;
                        velocities[counter] = Vector3.zero;

                        if (++counter == numberOfParticles)
                        {
                            return;
                        }
                    }
        }
    }

    private void InitNeighbourHashing()
    {
        _hashGrid.Clear();  // Only needed when resetting the simulation via testBattery approach.
        _dynamicCollisionHashGrid.Clear();
        _staticCollisionHashGrid.Clear();
        _neighbourList = new int[numberOfParticles * maximumParticlesPerCell * 8];   // 8 because we consider 8 cells
        _neighbourCollisionList = new Triangle[numberOfParticles * maximumParticlesPerCell * 8];   // 8 because we consider 8 cellsc
        _neighbourTracker = new int[numberOfParticles];
        _neighbourCollisionTracker = new int[numberOfParticles];
        HashGrid.CellSize = radius * 2; // Setting cell-size h to particle diameter.
        HashGrid.Dimensions = dimensions; 
        for (int i = 0; i < dimensions; i++)
            for (int j = 0; j < dimensions; j++)
                for (int k = 0; k < dimensions; k++)
                {
                    _hashGrid.Add(HashGrid.Hash(new Vector3Int(i, j, k)), new List<int>());
                    _dynamicCollisionHashGrid.Add(HashGrid.Hash(new Vector3Int(i, j, k)), new HashSet<Triangle>());
                    _staticCollisionHashGrid.Add(HashGrid.Hash(new Vector3Int(i, j, k)), new HashSet<Triangle>());
                }
    }

    #endregion


    // for calculations
    Vector3[] vertices;
    int[] triangles;
    Vector3[] normals;
    Matrix4x4 localToWorld;

    void Update()
    {
        // Calculate hash of all particles and build neighboring list.
        // 1. Clear HashGrid
        foreach (var cell in _hashGrid)
        {
            cell.Value.Clear();
        }
        foreach (var cell in _dynamicCollisionHashGrid)
        {
            cell.Value.Clear();
        }
        // 2. Recalculate hashes of each particle.
        for (int i = 0; i < _particles.Length; i++)
        {
            var hash = HashGrid.Hash(HashGrid.GetCell(_particles[i].transform.position));
            if (_hashGrid[hash].Count == maximumParticlesPerCell) continue;   // Prevent potential UB in neighbourList if more than maxParticlesPerCell are in a cell.
            _hashGrid[hash].Add(i);
        }
        // 2. Recalculate hashes of each collider
        InitColliders(dynamicColliders, false);

        // 3. For each particle go through all their 8 neighbouring cells.
        //    Check each particle in those neighbouring cells for interference radius r and store the interfering ones inside the particles neighbour list.
        for (int particleIndex = 0; particleIndex < _particles.Length; particleIndex++)
        {
            _neighbourTracker[particleIndex] = 0;
            _neighbourCollisionTracker[particleIndex] = 0;
            var cell = HashGrid.GetCell(_particles[particleIndex].transform.position);
            var cells = GetNearbyKeys(cell, _particles[particleIndex].transform.position);

            for (int j = 0; j < cells.Length; j++)
            {
                if (_hashGrid.ContainsKey(cells[j]))
                {
                    var neighbourCell = _hashGrid[cells[j]];
                    foreach (var potentialNeighbour in neighbourCell)
                    {
                        if (potentialNeighbour == particleIndex) continue;

                        if ((_particles[potentialNeighbour].transform.position - _particles[particleIndex].transform.position).sqrMagnitude < radius2)
                        {
                            _neighbourList[particleIndex * maximumParticlesPerCell * 8 + _neighbourTracker[particleIndex]++] = potentialNeighbour;
                        }
                    }
                }
                if (_dynamicCollisionHashGrid.ContainsKey(cells[j]))
                {
                    var neighbourCell = _dynamicCollisionHashGrid[cells[j]];
                    foreach (var potentialNeighbour in neighbourCell)
                    {
                        _neighbourCollisionList[particleIndex * maximumParticlesPerCell * 8 + _neighbourCollisionTracker[particleIndex]++] = potentialNeighbour;
                    }
                }
                if (_staticCollisionHashGrid.ContainsKey(cells[j]))
                {
                    var neighbourCell = _staticCollisionHashGrid[cells[j]];
                    foreach (var potentialNeighbour in neighbourCell)
                    {
                        _neighbourCollisionList[particleIndex * maximumParticlesPerCell * 8 + _neighbourCollisionTracker[particleIndex]++] = potentialNeighbour;
                    }
                }
            }
        }

        ComputeDensityPressure();
        ComputeForces();
        Integrate();
        averageFPS = 1 / Time.deltaTime;
    }

    private void MakeCollision(int i)
    {
        Vector3 newPos = _particles[i].transform.position;
        for (int j = 0; j < _neighbourCollisionTracker[i]; j++)
        {
            Triangle face = _neighbourCollisionList[i * maximumParticlesPerCell * 8 + j];
            Vector3 closestPoint = face.ClosestPointToPlane(newPos);
            if (Vector3.Distance(newPos, closestPoint) > radius/2 || !face.PointInTriangle(closestPoint) || Vector3.Dot(face.normal, velocities[i].normalized) >=0)
                continue;
            float collisionDamping = 1.0f;
            Vector3 force = Vector3.Project(velocities[i], (closestPoint - newPos).normalized);

            if (force.magnitude > 1) collisionDamping *= 0.25f; 

            velocities[i] = Vector3.Reflect(velocities[i], face.normal.normalized) * collisionDamping ;
        }
    }
    private void Integrate()
    {
        for (int i = 0; i < numberOfParticles; i++)
        {
            // forward Euler integration
            
            Vector3 acceleration = forces[i] / densities[i];
            velocities[i] += acceleration * 0.01f;
            Vector3 newPos = _particles[i].transform.position;
            MakeCollision(i);
            //velocities[i] += Time.deltaTime * forces[i] / mass;

            newPos += 0.01f * velocities[i];

            // enforce boundary conditions
            if (newPos.x - float.Epsilon < 0.0f)
            {
                velocities[i].x *= damping;
                newPos.x = float.Epsilon;
            }
            else if (newPos.x + float.Epsilon > dimensions - 1f)
            {
                velocities[i].x *= damping;
                newPos.x = dimensions - 1 - float.Epsilon;
            }

            if (newPos.y - float.Epsilon < 0.0f)
            {
                velocities[i].y *= damping;
                newPos.y = float.Epsilon;
            }
            else if (newPos.y + float.Epsilon > dimensions - 1f)
            {
                velocities[i].y *= damping;
                newPos.y = dimensions - 1 - float.Epsilon;
            }

            if (newPos.z - float.Epsilon < 0.0f)
            {
                velocities[i].z *= damping;
                newPos.z = float.Epsilon;
            }
            else if (newPos.z + float.Epsilon > dimensions - 1f)
            {
                velocities[i].z *= damping;
                newPos.z = dimensions - 1 - float.Epsilon;
            }
            _particles[i].transform.position = newPos;
        }
    }

    private int[] GetNearbyKeys(Vector3Int originIndex, Vector3 position)
    {
        Vector3Int[] nearbyBucketIndices = new Vector3Int[8];
        for (int i = 0; i < 8; i++)
        {
            nearbyBucketIndices[i] = originIndex;
        }

        if ((originIndex.x + 0.5f) * HashGrid.CellSize <= position.x)
        {
            nearbyBucketIndices[4].x += 1;
            nearbyBucketIndices[5].x += 1;
            nearbyBucketIndices[6].x += 1;
            nearbyBucketIndices[7].x += 1;
        }
        else
        {
            nearbyBucketIndices[4].x -= 1;
            nearbyBucketIndices[5].x -= 1;
            nearbyBucketIndices[6].x -= 1;
            nearbyBucketIndices[7].x -= 1;
        }

        if ((originIndex.y + 0.5f) * HashGrid.CellSize <= position.y)
        {
            nearbyBucketIndices[2].y += 1;
            nearbyBucketIndices[3].y += 1;
            nearbyBucketIndices[6].y += 1;
            nearbyBucketIndices[7].y += 1;
        }
        else
        {
            nearbyBucketIndices[2].y -= 1;
            nearbyBucketIndices[3].y -= 1;
            nearbyBucketIndices[6].y -= 1;
            nearbyBucketIndices[7].y -= 1;
        }

        if ((originIndex.z + 0.5f) * HashGrid.CellSize <= position.z)
        {
            nearbyBucketIndices[1].z += 1;
            nearbyBucketIndices[3].z += 1;
            nearbyBucketIndices[5].z += 1;
            nearbyBucketIndices[7].z += 1;
        }
        else
        {
            nearbyBucketIndices[1].z -= 1;
            nearbyBucketIndices[3].z -= 1;
            nearbyBucketIndices[5].z -= 1;
            nearbyBucketIndices[7].z -= 1;
        }

        int[] nearbyKeys = new int[8];
        for (int i = 0; i < 8; i++)
        {
            nearbyKeys[i] = HashGrid.Hash(nearbyBucketIndices[i]);
        }

        return nearbyKeys;
    }

    void ComputeDensityPressure()
    {
        for (int i = 0; i < _particles.Length; i++)
        {
            Vector3 origin = _particles[i].transform.position;
            float sum = 0f;
            for (int j = 0; j < _neighbourTracker[i]; j++)
            {
                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];
                float distanceSqr = (origin - _particles[neighbourIndex].transform.position).sqrMagnitude;

                sum += k1 * Mathf.Pow((radius2 - distanceSqr),3); 
            }

            sum += k1 * Mathf.Pow(radius, 6);

            densities[i] = sum;

            // Compute pressure based on density
            pressures[i] = gasConstant * (densities[i] - restDensity);
        }
    }

    public void ComputeForces()
    {
        for (int i = 0; i < _particles.Length; i++)
        {
            forces[i] = Vector3.zero;
            for (int j = 0; j < _neighbourTracker[i]; j++)
            {
                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];

                float distance = (_particles[i].transform.position - _particles[neighbourIndex].transform.position).magnitude;
                Vector3 direction = (_particles[i].transform.position - _particles[neighbourIndex].transform.position).normalized;

                Vector3 pressureForce = k2 * -(direction) * (pressures[i] + pressures[neighbourIndex]) / (2 * densities[neighbourIndex]) * Mathf.Pow((radius - distance),2); 
                forces[i] += pressureForce;

                // Viscosity 
                Vector3 velocityDiff = velocities[neighbourIndex] - velocities[i];
                Vector3 viscoForce = k3 * (velocityDiff / densities[neighbourIndex]) * (radius - distance);
                forces[i] += viscoForce;
            }
            // Gravity
            /*            Vector3 center = new Vector3(dimensions / 2, dimensions / 2, dimensions / 2);
                        forces[i] += (center - _particles[i].transform.position).normalized * 9.81f ;*/
            forces[i] += g;
            if (tsunamiMode)
                forces[i] += TsunamiForces(_particles[i].transform.position);
        }
    }

    public Vector3 TsunamiForces (Vector3 position)
    {
        float halfD = dimensions / 2.0f;
        int quarter = 0;


        if (position.x > halfD) quarter++;
        if (position.y > halfD)
        {
            quarter++;
            if (position.x <= halfD)
                quarter++; 
        }

        float x_percentage = position.x  / dimensions;
        float y_percentage = position.y / dimensions;

        switch (quarter)
        {
            case 0:
                return new Vector3( 20, 0, 0); 
            case 1:
                return new Vector3( -5 *  x_percentage , 25 * (1 - y_percentage), 0);
            case 2:
                return new Vector3( -50 * x_percentage, -100  * y_percentage, 0);
            case 3:
                return new Vector3( 300 , -500 * y_percentage, 0);
        }

        return Vector3.zero;
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = Color.red;
        Gizmos.DrawWireCube(new Vector3(dimensions / 2f, dimensions / 2f, dimensions / 2f), Vector3.one * dimensions);
    }

}

