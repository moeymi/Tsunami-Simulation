using UnityEngine;
public class Triangle
{
    public readonly Plane plane;
    public readonly Vector3 p0 = Vector3.zero, p1 = Vector3.zero, p2;
    public readonly Vector3 normal = Vector3.zero;
    public readonly Vector3 center = Vector3.zero;
    public readonly float triArea;
    public bool isStatic;

    public Triangle(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 normal)
    {
        this.p0 = p0;
        this.p1 = p1;
        this.p2 = p2;
        this.normal = normal;
        triArea = Vector3.Cross(p2 - p0, p2 - p1).magnitude / 2;
        plane = new Plane(normal, p0);
    }
    public override bool Equals(object obj)
    {
        Triangle other = obj as Triangle;
        return other.GetHashCode() == other.GetHashCode();
    }
    public override int GetHashCode()
    {
        Vector3 c = (p1 + p2 + p2) / 3;
        Vector3Int center = new Vector3Int((int)c.x % 7, (int)c.y % 7, (int)c.z % 7);
        string h = center.x.ToString() + center.y.ToString() + center.z.ToString();
        return int.Parse(h);
    }

    public Vector3 ClosestPointToPlane(Vector3 p)
    {
        return plane.ClosestPointOnPlane(p);
    }

    public bool PointInTriangle(Vector3 p)
    {
        Vector3 d1;
        Vector3 d2;
        float sumArea = 0;
        d1 = p - p0;
        d2 = p - p1;
        sumArea += Vector3.Cross(d1, d2).magnitude / 2;

        d1 = p - p1;
        d2 = p - p2;
        sumArea += Vector3.Cross(d1, d2).magnitude / 2;

        d1 = p - p0;
        d2 = p - p2;
        sumArea += Vector3.Cross(d1, d2).magnitude / 2;

        if (Mathf.Abs(sumArea - triArea) < 0.2)
        {
            return true;
        }
        return false;
    }
}