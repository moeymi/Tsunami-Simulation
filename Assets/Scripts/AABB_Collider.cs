using UnityEditor;
using UnityEngine;

public class AABB_Collider : MonoBehaviour
{
    [SerializeField]
    AABB_Collider another;
    public Vector3 size
    {
        get; private set;
    }

    private void Awake()
    {
        size = transform.lossyScale/2;
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = Color.yellow;
        Matrix4x4 rotationMatrix = Matrix4x4.TRS(transform.position, transform.rotation, transform.lossyScale);
        Gizmos.matrix = rotationMatrix;

        Gizmos.DrawWireCube(Vector3.zero,  Vector3.one);
    }

    private void Update()
    {
        // collision x-axis?
        bool collisionX = another.transform.position.x + another.size.x >= transform.position.x - size.x &&
            transform.position.x + size.x >= another.transform.position.x - another.size.x;
        // collision y-axis?
        bool collisionY = another.transform.position.y + another.size.y >= transform.position.y - size.y &&
            transform.position.y + size.y >= another.transform.position.y - another.size.y;

        bool collisionZ = another.transform.position.z + another.size.z >= transform.position.z - size.z &&
            transform.position.z + size.z >= another.transform.position.z - another.size.z;
        if (collisionX && collisionY && collisionZ)
        {
            GetComponent<Renderer>().material.color = Color.red;
        }
        else
        {
            GetComponent<Renderer>().material.color = Color.white;
        }
    }
}
