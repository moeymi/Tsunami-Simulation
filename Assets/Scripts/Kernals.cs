using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Kernels
{
	static public float Poly6(float distSqr, float h)
	{
		float coef = 315f / (64f * Mathf.PI * Mathf.Pow(h, 9));
		float hSqr = h * h;

		if (hSqr < distSqr)
			return 0f;

		return coef * Mathf.Pow(hSqr - distSqr, 3);
	}

	static public Vector2 GradientSpiky(Vector2 r, float h)
	{
		float coef = 45f / (Mathf.PI * Mathf.Pow(h, 6));
		float dist = r.magnitude;

		if (h < dist)
			return Vector3.zero;

		return -coef * r.normalized * Mathf.Pow(h - dist, 2);
	}

	static public float ViscosityLaplacian(float r, float h)
	{
		if (h < r)
			return 0f;

		float coef = 45f / (Mathf.PI * Mathf.Pow(h, 6));
		return coef * (h - r);
	}


	static public float Poly6New(float h)
	{
		float coef = 315f / (64f * Mathf.PI * Mathf.Pow(h, 9));
		float h2 = h * h;  
		return coef ; 
	}

	static public float GradientSpikyNew(float r,float h)
	{
		if (h < r)
			return 0f;
		float coef = 15f / (Mathf.PI * Mathf.Pow(h, 6));
		return coef * Mathf.Pow(h  - r , 3); 
	}
	static public float ViscosityLaplacianNew(float r, float h)
	{
		if (h < r)
			return 0f;

		float coef = 45f / (Mathf.PI * Mathf.Pow(h, 6));
		return coef * ( h - r );
	}
}
