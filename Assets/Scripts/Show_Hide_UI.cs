using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Show_Hide_UI : MonoBehaviour
{
    public GameObject panel1, panel2; 
    public void show_hide ()
    {
        bool state = !panel1.activeSelf;
        panel1.SetActive(state);
        panel2.SetActive(state);
    }
}
