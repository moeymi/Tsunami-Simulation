using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEditor;
using UnityEngine;
using Debug = UnityEngine.Debug;
using Random = UnityEngine.Random;

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
    public int numberOfParticles = 1000;
    public int dimensions = 10;
    public int maximumParticlesPerCell = 500;

    [Header("Debug information")]
    [Tooltip("Tracks how many neighbours each particleIndex has in " + nameof(_neighbourList))]
    [SerializeField]
    private int[] _neighbourTracker;

    public float averageFPS;
    #endregion

    private GameObject[] _particles;
    private int[] _neighbourList; // Stores all neighbours of a particle aligned at 'particleIndex * maximumParticlesPerCell * 8'
    private readonly Dictionary<int, List<int>> _hashGrid = new Dictionary<int, List<int>>();  // Hash of cell to particle indices.

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

    float minVal = float.MaxValue;
    float maxVal = float.MinValue;

    private void Awake()
    {
        RespawnParticles();
        InitNeighbourHashing();
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;

        k1 = (315 * mass) / (64 * Mathf.PI * Mathf.Pow(radius, 9));
        k2 = -(45 * mass) / (Mathf.PI * Mathf.Pow(radius, 6));
        k3 = (45 * viscosityCoefficient * mass) / (Mathf.PI * Mathf.Pow(radius, 6)); 
    }

    #region Initialisation

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
        while (counter < numberOfParticles)
        {
            for (int x = 0; x < particlesPerDimension; x++)
                for (int y = 0; y < particlesPerDimension; y++)
                    for (int z = 0; z < particlesPerDimension; z++)
                    {
                        //Vector3 startPos = new Vector3(dimensions - 1, dimensions - 1, dimensions - 1) - new Vector3(x / 2f, y / 2f, z / 2f) - new Vector3(Random.Range(0f, 1f), Random.Range(0f, 1f), Random.Range(0f, 1f));
                        Vector3 startPos = new Vector3(dimensions -1 - x, dimensions -1 - y, dimensions -1 - z); 
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
        _neighbourList = new int[numberOfParticles * maximumParticlesPerCell * 8];   // 8 because we consider 8 cells
        _neighbourTracker = new int[numberOfParticles];
        HashGrid.CellSize = radius * 2; // Setting cell-size h to particle diameter.
        HashGrid.Dimensions = dimensions; 
        for (int i = 0; i < dimensions; i++)
            for (int j = 0; j < dimensions; j++)
                for (int k = 0; k < dimensions; k++)
                {
                    _hashGrid.Add(HashGrid.Hash(new Vector3Int(i, j, k)), new List<int>());
                }
    }

    #endregion

    void Update()
    {
        // Calculate hash of all particles and build neighboring list.
        // 1. Clear HashGrid
        foreach (var cell in _hashGrid)
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
        // 3. For each particle go through all their 8 neighbouring cells.
        //    Check each particle in those neighbouring cells for interference radius r and store the interfering ones inside the particles neighbour list.
        for (int particleIndex = 0; particleIndex < _particles.Length; particleIndex++)
        {
            _neighbourTracker[particleIndex] = 0;
            var cell = HashGrid.GetCell(_particles[particleIndex].transform.position);
            var cells = GetNearbyKeys(cell, _particles[particleIndex].transform.position);

            for (int j = 0; j < cells.Length; j++)
            {
                if (!_hashGrid.ContainsKey(cells[j])) continue;
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
        }
        // 4. The Neighbouring-list should be n-particles big, each index containing a list of each particles neighbours in radius r.

        ComputeDensityPressure();
        ComputeForces();
        Integrate();
        averageFPS = 1 / Time.deltaTime;
    }

    private void Integrate()
    {
        for (int i = 0; i < numberOfParticles; i++)
        {
            // forward Euler integration
            
            Vector3 acceleration = forces[i] / densities[i];
            velocities[i] += acceleration *0.01f;
            //velocities[i] += Time.deltaTime * forces[i] / mass;

            Vector3 newPos = _particles[i].transform.position;
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
            forces[i] += g;
        }
    }

    #region OLD COMPUTE METHODS
    private void ComputeForces2()
    {
        float mass2 = mass * mass;
        Vector3 CS_D1 = Vector3.zero;
        float CS_D2 = 0f;
        float K = 0f;
        float surfaceTensionCo = (72 / 300f);

        for (int i = 0; i < _particles.Length; i++)
        {
            forces[i] = Vector3.zero;
            var particleDensity2 = densities[i] * densities[i];
            Vector3 viscosityForce = Vector3.zero;
            Vector3 surfaceTensionForce = Vector3.zero;
            for (int j = 0; j < _neighbourTracker[i]; j++)
            {
                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];
                float distance = (_particles[i].transform.position - _particles[neighbourIndex].transform.position).magnitude;
                Debug.Log(i + " distance = " + distance);
                if (distance > 0)
                {
                    var direction = (_particles[i].transform.position - _particles[neighbourIndex].transform.position).normalized;

                    //Vector3 kernel = SpikyKernelGradient(distance, direction);
                    Vector3 kernel = Kernels.GradientSpiky((_particles[i].transform.position - _particles[neighbourIndex].transform.position), radius);
                    float bigVal = (pressures[i] / particleDensity2 + pressures[neighbourIndex] / (densities[neighbourIndex] * densities[neighbourIndex]));
                    forces[i] -= mass2 * bigVal * kernel;
                    Debug.Log(i + " Density = " + densities[i]);
                    Debug.Log(i + " Pressures = " + forces[i]);

                    // 8. Compute the viscosity force
                    Vector3 visc = (velocities[neighbourIndex] - velocities[i]) / densities[neighbourIndex];
                    //float kernelVis = SpikyKernelSecondDerivative(distance);
                    float kernelVis = Kernels.ViscosityLaplacian(distance, radius);
                    viscosityForce += mass * visc * kernelVis;  // Kim

                    //Compute Surface Tension 
                    CS_D1 = mass * (1 / densities[neighbourIndex]) * Kernels.GradientSpiky((_particles[i].transform.position - _particles[neighbourIndex].transform.position), radius);
                    CS_D2 = mass * (1 / densities[neighbourIndex]) * Kernels.ViscosityLaplacian(distance, radius);
                    K = -CS_D2 / CS_D1.magnitude;
                    if (CS_D1.magnitude > 0.2)
                        surfaceTensionForce += surfaceTensionCo * K * CS_D1;


                }
            }
            viscosityForce *= viscosityCoefficient;
            forces[i] += viscosityForce;
            forces[i] += surfaceTensionForce;

            // Gravity
            forces[i] += g;
        }
    }

    private void ComputeDensityPressure2()
    {
        for (int i = 0; i < _particles.Length; i++)
        {
            Vector3 origin = _particles[i].transform.position;
            float sum = 0f;
            for (int j = 0; j < _neighbourTracker[i]; j++)
            {
                /*                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];
                                float distance = (origin - _particles[neighbourIndex].transform.position).sqrMagnitude;*/
                //sum += mass * StdKernel(distance);
                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];
                float distanceSqr = (origin - _particles[neighbourIndex].transform.position).sqrMagnitude;
                sum += mass * Kernels.Poly6New(radius) * Mathf.Pow(distanceSqr - radius2, 2);
            }

            //self density
            sum += mass * Kernels.Poly6New(radius) * Mathf.Pow(radius, 6);

            densities[i] = sum;

            // 6. Compute pressure based on density
            pressures[i] = gasConstant * (densities[i] - restDensity);
        }
    }


    private void ComputeDensityPressure3()
    {
        for (int i = 0; i < _particles.Length; i++)
        {
            Vector3 origin = _particles[i].transform.position;
            float sum = 0f;
            for (int j = 0; j < _neighbourTracker[i]; j++)
            {
                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];
                float distanceSqr = (origin - _particles[neighbourIndex].transform.position).sqrMagnitude;
                sum += mass * Kernels.Poly6New(radius) * Mathf.Pow(distanceSqr - radius2, 2);
            }

            // self density
            sum += mass * Kernels.Poly6New(radius) * Mathf.Pow(radius, 6);

            densities[i] = sum;

            // Compute pressure based on density
            pressures[i] = gasConstant * (densities[i] - restDensity);
        }
    }

    public void ComputeForces3()
    {
        for (int i = 0; i < _particles.Length; i++)
        {
            forces[i] = Vector3.zero;
            for (int j = 0; j < _neighbourTracker[i]; j++)
            {
                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];
                float distance = (_particles[i].transform.position - _particles[neighbourIndex].transform.position).magnitude;
                Vector3 direction = (_particles[i].transform.position - _particles[neighbourIndex].transform.position).normalized;

                //Pressure force 
                Vector3 pressureForce = (-direction * mass) * (pressures[i] + pressures[neighbourIndex]) / (2 * densities[neighbourIndex]);
                pressureForce *= Kernels.GradientSpikyNew(distance, radius);
                forces[i] += pressureForce;

                // Viscosity 
                Vector3 velocityDiff = velocities[neighbourIndex] - velocities[i];
                Vector3 viscoForce = viscosityCoefficient * mass * (velocityDiff / densities[neighbourIndex]) * Kernels.ViscosityLaplacian(distance, radius);
                forces[i] += viscoForce;
            }
            // Gravity
            forces[i] += g;
        }
    }

    // Kernel by Müller et al.
    private float StdKernel(float distanceSquared)
    {
        // Doyub Kim
        float x = 1.0f - distanceSquared / radius2;
        return 315f / (64f * Mathf.PI * radius3) * x * x * x;
    }

    // Doyub Kim page 130
    private float SpikyKernelFirstDerivative(float distance)
    {
        float x = 1.0f - distance / radius;
        return -45.0f / (Mathf.PI * radius4) * x * x;
    }

    // Doyub Kim page 130
    private float SpikyKernelSecondDerivative(float distance)
    {
        // Btw, it derives 'distance' not 'radius' (h)
        float x = 1.0f - distance / radius;
        return 90f / (Mathf.PI * radius5) * x;
    }

    // Doyub Kim page 130
    private Vector3 SpikyKernelGradient(float distance, Vector3 directionFromCenter)
    {
        return SpikyKernelFirstDerivative(distance) * directionFromCenter;
    }

    #endregion



}

