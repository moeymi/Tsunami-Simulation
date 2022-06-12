using UnityEngine;

public class Show_Hide_UI : MonoBehaviour
{
    public GameObject panel1, panel2;

    [SerializeField]
    private Camera[] _cameras;

    public void show_hide ()
    {
        bool state = !panel1.activeSelf;
        panel1.SetActive(state);
        panel2.SetActive(state);
    }
    public void changeCamera()
    {
        _cameras[0].enabled = !_cameras[0].enabled;
        _cameras[1].enabled = !_cameras[1].enabled;
        RenderSettings.fog = !RenderSettings.fog;
    }
}
