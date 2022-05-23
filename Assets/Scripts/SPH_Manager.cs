using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEditor;
using UnityEngine;
using Debug = UnityEngine.Debug;
using Random = UnityEngine.Random;

public class SPH_Manager : MonoBehaviour
{
    // ReSharper disable InconsistentNaming
    [Header("Particle properties")]
    public float radius = 1f; // particle radius, interaction radius h
    public GameObject particlePrefab;
    public float particleRadius = 1f;
    public float mass = 4f;
    public float viscosityCoefficient = 2.5f;
    private static readonly Vector3 g = new Vector3(0.0f, -9.81f, 0.0f);
    private const float gasConstant = 8.315f;

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

    private GameObject[] _particles;
    private int[] _neighbourList; // Stores all neighbours of a particle aligned at 'particleIndex * maximumParticlesPerCell * 8'
    private readonly Dictionary<int, List<int>> _hashGrid = new Dictionary<int, List<int>>();  // Hash of cell to particle indices.

    private float radius2;
    private float radius3;
    private float radius4;
    private float radius5;
    // ReSharper restore InconsistentNaming

    [StructLayout(LayoutKind.Sequential, Size = 28)]
    private struct Particle
    {
        public Vector3 Position;
        public Vector4 Color;
    }

    private float[] densities;
    private float[] pressures;
    private Vector3[] forces;
    private Vector3[] velocities;

    [InitializeOnLoadMethod]
    static public void ggg()
    {

        double bigVal = -1270740000;
        Debug.Log(bigVal);
        Debug.Log(bigVal * 0.002483614);
    }

    private void Awake()
    {
        RespawnParticles();
        InitNeighbourHashing();
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
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
                        Vector3 startPos = new Vector3(dimensions - 1, dimensions - 1, dimensions - 1) - new Vector3(x / 2f, y / 2f, z / 2f) - new Vector3(Random.Range(0f, 1f), Random.Range(0f, 1f), Random.Range(0f, 1f));
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

                    if ((_particles[potentialNeighbour].transform.position - _particles[particleIndex].transform.position).magnitude < radius)
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

    // https://lucasschuermann.com/writing/implementing-sph-in-2d
    private void Integrate()
    {
        for (int i = 0; i < numberOfParticles; i++)
        {
            // forward Euler integration
            velocities[i] += Time.deltaTime * forces[i] / mass;
            Vector3 newPos = _particles[i].transform.position;
            newPos += Time.deltaTime * velocities[i];

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

    float minVal = float.MaxValue;
    float maxVal = float.MinValue;

    private void ComputeForces()
    {
        float mass2 = mass * mass;
        for (int i = 0; i < _particles.Length; i++)
        {
            forces[i] = Vector3.zero;
            var particleDensity2 = densities[i] * densities[i];
            Vector3 viscosityForce = Vector3.zero;
            Vector3 surfaceTension = Vector3.zero;
            for (int j = 0; j < _neighbourTracker[i]; j++)
            {
                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];
                float distance = (_particles[i].transform.position - _particles[neighbourIndex].transform.position).magnitude;
                if (distance > 0.0f)
                {
                    var direction = (_particles[i].transform.position - _particles[neighbourIndex].transform.position) / distance;
                    // 7. Compute pressure gradient force (Doyub Kim page 136)
                    Vector3 kernel = SpikyKernelGradient(distance, direction);
                    float bigVal = (pressures[i] / particleDensity2 + pressures[neighbourIndex] / (densities[neighbourIndex] * densities[neighbourIndex]));
                    forces[i] -= mass2 * bigVal * kernel;   // Kim

                    // 8. Compute the viscosity force
                    Vector3 visc = (velocities[neighbourIndex] - velocities[i]) / densities[neighbourIndex];
                    float kernelVis = SpikyKernelSecondDerivative(distance);
                    viscosityForce += mass * visc * kernelVis;  // Kim

                    surfaceTension += mass * (1 / densities[neighbourIndex]) * SpikyKernelGradient(distance, direction);

                }
            }
            viscosityForce *= viscosityCoefficient;
            forces[i] += viscosityForce;
            forces[i] += surfaceTension;

            // Gravity
            forces[i] += g;
            if (forces[i].magnitude != 0)
            {
                minVal = Mathf.Min(minVal, forces[i].magnitude);
                maxVal = Mathf.Max(maxVal, forces[i].magnitude);
            }
        }
        Debug.Log(minVal);
        Debug.Log(maxVal);
    }

    private void ComputeDensityPressure()
    {
        for (int i = 0; i < _particles.Length; i++)
        {
            // Doyub Kim 121, 122, 123
            // 5. Compute densities
            Vector3 origin = _particles[i].transform.position;
            float sum = 0f;
            for (int j = 0; j < _neighbourTracker[i]; j++)
            {
                int neighbourIndex = _neighbourList[i * maximumParticlesPerCell * 8 + j];
                float distance = (origin - _particles[neighbourIndex].transform.position).magnitude;
                sum += mass * StdKernel(distance);
            }
            if(sum == 0)
                densities[i] = 1;
            else
                densities[i] = sum + 0.000001f;

            // 6. Compute pressure based on density
            pressures[i] = gasConstant * (densities[i] - restDensity); // as described in Müller et al Equation 12
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

    // Derived from Doyub Kim
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
}

